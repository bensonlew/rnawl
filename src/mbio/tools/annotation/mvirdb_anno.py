# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.files.align.blast.blast_table import BlastTableFile
import subprocess
import os


class MvirdbAnnoAgent(Agent):
    """
    宏基因组比对结果文件注释,先将xml文件转换成table并合并成一张table，进行anno详细注释
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(MvirdbAnnoAgent, self).__init__(parent)
        options = [
            {"name": "mvirdb_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # 比对到mvirdb库的xml文件夹
            {"name": "mvirdb_anno_result", "type": "outfile", 'format': "sequence.profile_table"}  # 注释详细结果表
        ]
        self.add_option(options)
        self.step.add_steps("merge", "mvirdb_anno")
        self.result_name = ''
        # self.on('start', self.stepstart)
        # self.on('end', self.stepfinish)

    def check_options(self):
        if not self.option("mvirdb_xml_dir").is_set:
            raise OptionError("必须设置输入文件", code="31204901")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        self.option('mvirdb_anno_result', os.path.join(self.output_dir, "gene_mvirdb_anno.xls"))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['gene_mvirdb_anno.xls', 'xls', '序列详细注释文件']
        ])
        super(MvirdbAnnoAgent, self).end()


class MvirdbAnnoTool(Tool):
    def __init__(self, config):
        super(MvirdbAnnoTool, self).__init__(config)
        self.logger.info(">>>init super success")
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        # self.python_path = "program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/annotation/dna_mvirdb.py'
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.result_name = ''


    def run(self):
        """
        运行
        :return:
        """
        super(MvirdbAnnoTool, self).run()
        self.set_output()
        self.merge_table()
        self.run_mvirdb_anno()
        self.end()


    def merge_table(self):
        self.mvirdb_number = 0
        xml_file = os.listdir(self.option('mvirdb_xml_dir').prop['path'])
        self.result_name = os.path.join(self.output_dir, "mvirdb_align_table.xls")
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        for i in xml_file:
            self.mvirdb_number += 1
            file_path = os.path.join(self.option('mvirdb_xml_dir').prop['path'], i)
            table = xml2table(file_path,
                              self.work_dir + "/tmp_mvirdb_anno/" + "mvirdb_" + str(self.mvirdb_number) + "_table.xls")
            if self.mvirdb_number > 1:
                os.system("sed -i " + '1d ' + table)
            cmd = '{} {} {}'.format(self.sh_path, table, self.result_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat" + str(self.mvirdb_number)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i), code="31204901")
                self.set_error("cat %s error", variables=(i), code="31204902")

    def run_mvirdb_anno(self):
        b = BlastTableFile()
        b.set_path(self.result_name)
        new_result = self.output_dir + "/top1_mvirdb_align_table.xls"
        b.retain_highest_score(new_result)
        cmd = '{} {} -i {} -o {}'.format(self.python_path, self.python_script, new_result,
                                         self.output_dir + "/gene_mvirdb_anno.xls")
        # command = self.add_command("anno", cmd).run()
        self.logger.info(cmd)
        # self.wait(command)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行mvirdb_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行mvirdb_anno出错', code="31204903")

    def set_output(self):
        if os.path.exists(self.work_dir + '/tmp_mvirdb_anno'):
            pass
        else:
            os.mkdir(self.work_dir + '/tmp_mvirdb_anno')
            # self.option('mvirdb_anno_result',os.path.join(self.output_dir,"gene_mvirdb_anno.xls"))
            # if len(os.listdir(self.output_dir)) == 1:
            #    self.logger.info("output right")
