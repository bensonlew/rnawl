# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import os


class PhiAnnoAgent(Agent):
    """
    宏基因组比对结果文件注释,先将xml文件转换成table并合并成一张table，进行anno详细注释
    author: haidong.gu
    last_modify:
    """

    def __init__(self, parent):
        super(PhiAnnoAgent, self).__init__(parent)
        options = [
            {"name": "phi_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # 比对到phi库的xml文件夹
            {"name": "phi_anno_result", "type": "outfile", 'format': "sequence.profile_table"},  # 注释详细结果表
            {"name": "project", "type": "string", "default": "dna"}  # 基因组还是宏基因组  dna/meta
            ]
        self.add_option(options)
        #self.phi_number = 0
        self.step.add_steps("merge","phi_anno")
        self.result_name = ''
        #self.on('start', self.stepstart)
        #self.on('end', self.stepfinish)

    def check_options(self):
        if not self.option("phi_xml_dir").is_set:
            raise OptionError("必须设置输入文件", code="31203301")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        self.option('phi_anno_result',os.path.join(self.output_dir,"gene_phi_anno.xls"))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['gene_phi_anno.xls', 'xls', '序列详细注释文件']
            ])
        super(PhiAnnoAgent, self).end()


class PhiAnnoTool(Tool):
    def __init__(self, config):
        super(PhiAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        #self.python_path = "program/Python/bin/python"
        try:
            if self.option("project") == "meta":
                self.python_script = self.config.PACKAGE_DIR + '/annotation/meta_phi.py'
            elif self.option("project") == "dna":
                self.python_script = self.config.PACKAGE_DIR + '/annotation/dna_phi.py'
        except:
            self.python_script = self.config.PACKAGE_DIR + '/annotation/dna_phi.py'
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.result_name = ''

    def run(self):
        """
        运行
        :return:
        """
        super(PhiAnnoTool, self).run()
        self.make_dir()
        self.merge_table()
        self.run_phi_anno()
        self.end()

    def merge_table(self):
        self.phi_number = 0
        xml_file = os.listdir(self.option('phi_xml_dir').prop['path'])
        self.result_name = os.path.join(self.output_dir, "phi_align_table.xls")
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        for i in xml_file:
            self.phi_number += 1
            file_path = os.path.join(self.option('phi_xml_dir').prop['path'], i)
            table = xml2table(file_path, self.work_dir + "/tmp_phi_anno/" + "phi_" + str(self.phi_number) + "_table.xls")
            if self.phi_number > 1:
                            os.system("sed -i " +  '1d ' + table)
            cmd = '{} {} {}'.format(self.sh_path, table, self.result_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat" + str(self.phi_number)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i), code="31203301")
                self.set_error("cat %s error", variables=(i), code="31203302")

    def run_phi_anno(self):
        cmd = '{} {} -i {} -o {}'.format(self.python_path, self.python_script,self.result_name , self.output_dir + "/gene_phi_anno.xls")
        #command = self.add_command("anno", cmd).run()
        self.logger.info(cmd)
        #self.wait(command)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行phi_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行phi_anno出错', code="31203303")

    def make_dir(self):
        if os.path.exists(self.work_dir + '/tmp_phi_anno'):
                pass
        else:
                os.mkdir(self.work_dir + '/tmp_phi_anno')
        #self.option('phi_anno_result',os.path.join(self.output_dir,"gene_phi_anno.xls"))
        #if len(os.listdir(self.output_dir)) == 1:
        #    self.logger.info("output right")


