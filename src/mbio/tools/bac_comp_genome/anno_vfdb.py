# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
#20190917

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import os


class AnnoVfdbAgent(Agent):
    """
    宏基因组vfdb数据库比对结果文件注释,先将xml文件转换成table并合并成一张table，进行anno详细注释
    """

    def __init__(self, parent):
        super(AnnoVfdbAgent, self).__init__(parent)
        options = [
            {"name": "vfdb_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # 比对到vfdb库的xml文件夹
            {"name": "database", "type": "string", "default": "core"},               ### 注释核心库还是预测库
            {"name": "vfdb_anno_result", "type": "outfile", 'format': "sequence.profile_table"},  # 注释详细结果表
            ]
        self.add_option(options)
        self.result_name = ''

    def check_options(self):
        if not self.option("vfdb_xml_dir").is_set:
            raise OptionError("必须设置输入文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AnnoVfdbAgent, self).end()


class AnnoVfdbTool(Tool):
    def __init__(self, config):
        super(AnnoVfdbTool, self).__init__(config)
        self._version = "1.0"
        self.python_path ="/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/bac_comp_genome/anno_vfdb.py'
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.ref = self.config.SOFTWARE_DIR + "/database/VFDB/vfdb_all.xls"
        self.result_name = ''

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoVfdbTool, self).run()
        self.merge_table()
        self.run_vfdb_anno()
        self.end()

    def merge_table(self):
        if os.path.exists(self.work_dir + '/tmp_vfdb_anno'):
            pass
        else:
            os.mkdir(self.work_dir + '/tmp_vfdb_anno')
        self.vfdb_number = 0
        xml_file = os.listdir(self.option('vfdb_xml_dir').prop['path'])
        vfdb_align_table = "vfdb_" + self.option("database") + "_align_table.xls"
        self.result_name = os.path.join(self.output_dir, vfdb_align_table)
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        for i in xml_file:
            self.vfdb_number += 1
            file_path = os.path.join(self.option('vfdb_xml_dir').prop['path'], i)
            self.logger.info("转换表格")
            if self.option("database") == "core":
                self.logger.info(file_path)
                table = xml2table(file_path, self.work_dir + "/tmp_vfdb_anno/" + "vfdb_" + str(
                    self.vfdb_number) + "_core_table.xls")
                self.logger.info(file_path)
            elif self.option("database") == "predict":
                table = xml2table(file_path, self.work_dir + "/tmp_vfdb_anno/" + "vfdb_" + str(
                    self.vfdb_number) + "_predict_table.xls")
            cmd = '{} {} {}'.format(self.sh_path, table, self.result_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat" + str(self.vfdb_number)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i), code="31203701")
                raise Exception("cat {} error".format(i))
        with open(self.result_name, "r") as f:
            head = f.next().strip()
        os.system("sed -i '/^Score\t/d ' " + self.result_name)
        os.system("sed -i '1i" + head + "\' " + self.result_name)

    def run_vfdb_anno(self):
        cmd = ''
        if self.option("database") == "core":
            cmd = '{} {} -s {} -i {} -o {}'.format(self.python_path, self.python_script, self.ref, self.result_name , self.output_dir + "/gene_vfdb_core_anno.xls")
        elif self.option("database") == "predict":
            cmd = '{} {} -s {} -i {} -o {}'.format(self.python_path, self.python_script, self.ref, self.result_name , self.output_dir + "/gene_vfdb_predict_anno.xls")
        command = self.add_command("anno", cmd).run()
        self.logger.info(cmd)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('运行vfdb_anno完成')
        else:
            self.set_error('运行vfdb_anno出错')

    def set_output(self):
        pass