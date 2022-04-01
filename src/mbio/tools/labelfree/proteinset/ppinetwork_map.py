## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdongxuan"
#last_modify:20170418

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re


class PpinetworkMapAgent(Agent):
    """
    调用map.r脚本，进行将基因ID mapping 到STRINGid
    version v1.0
    author: hongdongxuan
    last_modify: 20170418
    """
    def __init__(self, parent):
        super(PpinetworkMapAgent, self).__init__(parent)
        options = [
            {"name": "diff_exp_gene", "type": "infile", "format": "labelfree.common"},
            {"name": "species", "type": "int", "default": 9606},
            {"name": "seq_mapped", "type": "outfile", "format": "labelfree.common"},
            {"name": "unmapped_seq", "type": "outfile", "format": "labelfree.common"},
            {"name": "unmapped_db", "type": "outfile", "format": "labelfree.common"},
        ]
        self.add_option(options)
        self.step.add_steps("map")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.map.start()
        self.step.update()

    def stepfinish(self):
        self.step.map.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        species_list = [30611, 9598, 61853, 9593, 9606, 9544, 9483, 30608, 9601, 9478, 10141, 10020, 10090, 9986, 10116,
                        43179, 37347, 9685, 9913, 9739, 9669, 9796, 132908, 59463, 9646, 9823, 9785, 9813, 9371, 9361,
                        28377, 9031, 13735, 9103, 59729, 8049, 31033, 8090, 8083, 69293, 99883, 8128, 7955, 13616, 9258,
                        9305, 9315, 7897, 7757, 7719, 51511, 6239, 7227, 4932, 15368, 4513, 4641, 4533, 4538, 4555,
                        4558, 4577, 59689, 3702, 3711, 3847, 3694, 4081, 4113, 29760, 88036, 3218, 3055, 45157]
        if not self.option("diff_exp_gene").is_set:
            raise OptionError("必须输入含有gene_id的差异基因列表！")
        # if int(self.option('species')) not in species_list:
        #    raise OptionError("不能进行蛋白质互作分析，因为string数据库中不存在该物种的蛋白质互作组数据！")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '8G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["diff_exp_mapped.txt", "txt", "含有STRINGid结果信息"],

        ])
        super(PpinetworkMapAgent, self).end()


class PpinetworkMapTool(Tool):
    """
    将基因ID mapping 到STRINGid tool
    """
    def __init__(self, config):
        super(PpinetworkMapTool, self).__init__(config)
        self._version = '1.0.1'
        self.r_path = 'program/R-3.3.1/bin/Rscript'
        self.database = self.config.SOFTWARE_DIR + '/database/Annotation/all/String/string11.5/'
        # self.database = self.config.SOFTWARE_DIR + '/database/stringDB_10.5/'
        self.script_path = os.path.join(Config().PACKAGE_DIR, "labelfree/new_map.r")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/lib')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/lib64')

    def run_map(self):
        one_cmd = self.r_path + " %s %s %s %s %s" % (self.script_path, self.option('diff_exp_gene').prop['path'],
                                                           self.option('species'), 'PPI_result', self.database)
        self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功")
        else:
            self.logger.info("运行one_cmd出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir + '/PPI_result/')
        for f in results:
            if re.search(r'.*Rdata$', f):
                pass
            else:
                os.link(self.work_dir + '/PPI_result/' + f, self.output_dir + '/' + f)
            if f == 'seq_mapped.txt':
                self.option('seq_mapped').set_path(self.output_dir + '/' + f)
            if f == 'unmapped_seq.txt':
                self.option('unmapped_seq').set_path(self.output_dir + '/' + f)
            if f == 'unmapped_db.txt':
                self.option('unmapped_db').set_path(self.output_dir + '/' + f)
        self.logger.info('设置文件夹路径成功')


    def run(self):
        super(PpinetworkMapTool, self).run()
        self.run_map()
        self.set_output()
        self.end()
