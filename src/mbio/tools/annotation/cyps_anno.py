# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import os


class CypsAnnoAgent(Agent):
    """
    宏基因组交互分析，对比对结果文件进行注释，
    首先先将xml文件转成table，并合成一张table，再进行注释
    """
    def __init__(self, parent):
        super(CypsAnnoAgent, self).__init__(parent)
        options = [
            {"name": "cyps_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # 比对到p450库的xml文件夹
            {"name": "cyps_anno_result", "type": "outfile", "format": "sequence.profile_table"}, #注释详细结果表
            {"name": "result_format", "type": "string", "default": "meta"}  #输出结果类型
        ]
        self.add_option(options)
        self.step.add_steps("merge", "cyps_anno")
        self.result_name = ''

    def check_options(self):
        if not self.option("cyps_xml_dir").is_set:
            raise OptionError("必须设置输入文件", code="31204701")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        self.option("cyps_anno_result", os.path.join(self.output_dir, "gene_cyps_anno.xls"))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["gene_cyps_anno.xls", "xls", "序列详细注释文件"]
        ])
        super(CypsAnnoAgent, self).end()


class CypsAnnoTool(Tool):
    def __init__(self, config):
        super(CypsAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/annotation/scripts/meta_cyps_mongo.py'
        #self.python_script ='/mnt/ilustre/users/sanger-dev/sg-users/zhangqingchen/test/cyps_test/package/meta_cyps_mongo.py'
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.result_name = ''

    def run(self):
        """
        运行tool
        :return:
        """
        super(CypsAnnoTool, self).run()
        self.set_output()
        self.merge_table()
        self.run_cyps_anno()
        self.end()

    def merge_table(self):
        self.cyps_number = 0
        xml_file = os.listdir(self.option('cyps_xml_dir').prop["path"])
        self.result_name = os.path.join(self.output_dir, "cyps_align_table.xls")
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        for i in xml_file:
            self.cyps_number += 1
            file_path = os.path.join(self.option('cyps_xml_dir').prop['path'], i)
            table = xml2table(file_path, self.work_dir + "/tmp_cyps_anno" + "cyps" + str(self.cyps_number) + "_table.xls")
            cmd = '{} {} {}'.format(self.sh_path, table, self.result_name)
            self.logger.info("start cat{}".format(i))
            command_name = "cat" + str(self.cyps_number)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i), code="31204701")
                self.set_error("cat %s error", variables=(i), code="31204702")
        with open(self.result_name, "r") as f:
            head = f.next().strip()
        os.system("sed -i '/^Score\t/d' " + self.result_name)
        os.system("sed -i '1i'" + head + "\'" + self.result_name)

    def run_cyps_anno(self):
        cmd = '{} {} -i {} -o {}'.format(self.python_path, self.python_script, self.result_name, self.output_dir + "/gene_cyps_anno.xls")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行cyps_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行cyps_anno出错', code="31204703")

    def set_output(self):
        if os.path.exists(self.work_dir + "/tmp_cyps_anno"):
            pass
        else:
            os.mkdir(self.work_dir + "/tmp_cyps_anno")
        self.logger.info('cyps_anno 注释完成，正在输出结果')





















