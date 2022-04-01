     # -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.files.align.blast.blast_table import BlastTableFile
import subprocess
import os


class AnnoTcdbAgent(Agent):
    """
    宏基因组比对结果文件注释,先将xml文件转换成table并合并成一张table，进行anno详细注释
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(AnnoTcdbAgent, self).__init__(parent)
        options = [
            {"name": "tcdb_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # 比对到tcdb库的xml文件夹
            {"name": "tcdb_anno_result", "type": "outfile", 'format': "sequence.profile_table"},  # 注释详细结果表
        ]
        self.add_option(options)
        self.step.add_steps("merge", "tcdb_anno")

    def check_options(self):
        if not self.option("tcdb_xml_dir").is_set:
            raise OptionError("必须设置输入文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AnnoTcdbAgent, self).end()


class AnnoTcdbTool(Tool):
    def __init__(self, config):
        super(AnnoTcdbTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        self.db = self.config.SOFTWARE_DIR + '/database/TCDB/TCDB_function.xls'
        self.python_script = self.config.PACKAGE_DIR + '/bac_comp_genome/anno_tcdb.py'
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.result_name = ''

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoTcdbTool, self).run()
        self.set_output()
        self.merge_table()
        self.run_tcdb_anno()
        self.end()

    def merge_table(self):
        self.tcdb_number = 0
        xml_file = os.listdir(self.option('tcdb_xml_dir').prop['path'])
        self.result_name = os.path.join(self.output_dir, "tcdb_align_table.xls")
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        for i in xml_file:
            self.tcdb_number += 1
            file_path = os.path.join(self.option('tcdb_xml_dir').prop['path'], i)
            table = xml2table(file_path,
                              self.work_dir + "/tmp_tcdb_anno/" + "tcdb_" + str(self.tcdb_number) + "_table.xls")
            if self.tcdb_number > 1:
                os.system("sed -i " + '1d ' + table)
            cmd = '{} {} {}'.format(self.sh_path, table, self.result_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat" + str(self.tcdb_number)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i))

    def run_tcdb_anno(self):
        b = BlastTableFile()
        b.set_path(self.result_name)
        new_result = self.output_dir + "/top1_tcdb_align_table.xls"
        b.retain_highest_score(new_result)
        cmd = '{} {} -f {} -i {} -o {}'.format(self.python_path, self.python_script, self.db,  new_result,
                                         self.output_dir + "/gene_tcdb_anno.xls")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行tcdb_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行tcdb_anno出错')

    def set_output(self):
        if os.path.exists(self.work_dir + '/tmp_tcdb_anno'):
            pass
        else:
            os.mkdir(self.work_dir + '/tmp_tcdb_anno')
