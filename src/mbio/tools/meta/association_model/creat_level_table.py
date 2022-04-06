# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import unittest

class CreatLevelTableAgent(Agent):
    """
    根据gene_anno和gene_profile计算某一层级所有丰度并取top
    author: shaohua.yuan
    last_modify: haidong.gu
    """

    def __init__(self, parent):
        super(CreatLevelTableAgent, self).__init__(parent)
        options = [
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_list", "type": "infile", 'format': "sequence.profile_table"},
            {"name": "anno_file", "type": "infile", 'format': "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "level", "type": "string"},
            {"name": "top", "type": "string", "default": "all"},
            {"name": "group_method", "type": "int", "default": 0},
            {"name": "outprofile", "type": "outfile", 'format': "sequence.profile_table"},
            {"name": "anno_type", "type": "string"},     # nr时使用
            {"name": "samples", "type": "string"},
            {"name": "Total", "type": "int", "default": 1},   ## creat的表格结果是否需要Total列,1:need or 0:remove
            {"name": "levelname", "type": "string"},
            {"name": "lowestlevel", "type": "string"},
            {"name": "colorlevel", "type": "string"},
            {"name": "database", "type": "string"},
        ]
        self.add_option(options)
        self.outfile = ''
        self._memory_increase_step = 50  # 每次重运行增加50G内存 add by GHD @ 20180711

    def check_options(self):
        if not self.option("gene_profile").is_set:
            raise OptionError("必须设置输入gene_profile文件", code="32701201")
        if not self.option("gene_list").is_set:
            raise OptionError("必须设置输入gene_list文件", code="32701202")
        if not self.option("anno_file").is_set:
            raise OptionError("必须设置输入注释文件", code="32701203")
        if not (self.option("Total") == 0 or self.option("Total") == 1):
            raise OptionError("Total参数必须为0或1", code="32701204")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'  # 改回 by guhaidong @ 20180427
        # memory = 5 + 10 * self._rerun_time  # 每次重运行增加5G内存 by guhaidong @ 20180417
        # self._memory = "%sG" % memory

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(CreatLevelTableAgent, self).end()


class CreatLevelTableTool(Tool):
    def __init__(self, config):
        super(CreatLevelTableTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/metagenomic/scripts/creat_level_table.py'

    def run(self):
        """
        运行
        :return:
        """
        super(CreatLevelTableTool, self).run()
        self.run_select()
        self.set_output()
        self.end()

    def run_select(self):
        anno_file = self.option("anno_file").prop["path"]
        gene_list = self.option("gene_list").prop["path"]
        gene_profile = self.option("gene_profile").prop["path"]
        self.outfile = os.path.join(self.output_dir, "select_level_abu.xls")
        cmd = "{} {} -i {} -g {} -p {} -l '{}' -t {} -gm {} -total {} -o {}".format(self.python_path, self.script,
                                                                        anno_file,
                                                                        gene_list, gene_profile,
                                                                        self.option("level"),
                                                                        self.option("top"),
                                                                        self.option("group_method"),self.option("Total"),
                                                                        self.outfile)
        if self.option("group_table").is_set:
            group_file = self.option("group_table").prop["path"]
            cmd += " -m " + group_file
        if self.option("anno_type") == "nr":
            cmd += " -database nr"
        if self.option("samples"):
            cmd += " -sam " + self.option("samples")
        if "levelname" in self.get_option_object().keys() and self.option("levelname"):
            cmd += " -ln '{}' ".format(self.option("levelname"))
        if "lowestlevel" in self.get_option_object().keys() and self.option("lowestlevel"):
            cmd += " -lowest '{}' ".format(self.option("lowestlevel"))
        if "colorlevel" in self.get_option_object().keys() and self.option("colorlevel"):
            cmd += " -cl '{}' ".format(self.option("colorlevel"))
        if "database" in self.get_option_object().keys() and self.option("database") and self.option("database") != "nr":
            cmd += " -database " + self.option("database")
        command = self.add_command("creat_level_table", cmd, ignore_error=True).run()
        self.logger.info(cmd)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("creat_level_table succeed")
        elif command.return_code in [1, -9, -11]:  # add memory limit by guhaidong @ 20180417
            self.check_cmd("creat_level_table.o")
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("creat_level_table failed", code="32701201")
            self.set_error("creat_level_table failed", code="32701202")

    def check_cmd(self, log_file):
        no_result = None
        with open(log_file, 'r') as r:
            for line in r:
                if "参数下数据为空，请重新设置该参数!" in line:
                    no_result = line
        if no_result:
            self.set_error(line)

    def set_output(self):
        try:
            self.option('outprofile', self.outfile)
        except Exception as e:
            self.set_error("设置输出新level丰度表失败——%s", variables=(e), code="32701203")

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "CreatLevelTable",
            "type": "tool",
            "name": "meta.association_model.creat_level_table",
            "instant": True,
            "options": dict(
                gene_profile="/mnt/ilustre/users/sanger-dev/workspace/20181210/AnnoGo_tsg_31796_24475_807043/remote_input/gene_profile/reads_number.xls",
                gene_list="/mnt/ilustre/users/sanger-dev/workspace/20181210/AnnoGo_tsg_31796_24475_807043/remote_input/geneset_table/gene.uniGeneset.fa.length.txt",
                anno_file="/mnt/ilustre/users/sanger-dev/workspace/20181210/AnnoGo_tsg_31796_24475_807043/MgFunSelect/AnnoFunSelect/output/function_select_anno.xls",
                #group_table="/mnt/ilustre/users/sanger-dev/workspace/20181115/Metabolome_tsg_32888/noQC_group.xls",
                level="GO Term (Lev4)",
                top="50",
                samples="N2,N1,N4,N3,S1,N5,S3,S2,S5,S4,C1,C2,C3,C4,C5"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

