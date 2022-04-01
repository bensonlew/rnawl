# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import os


class QsAnnoAgent(Agent):
    """
    宏基因组比对结果文件注释,先将xml文件转换成table并合并成一张table，进行anno详细注释
    last_modify:
    """

    def __init__(self, parent):
        super(QsAnnoAgent, self).__init__(parent)
        options = [
            {"name": "qs_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # 比对到QS库的xml文件夹
            {"name": "qs_anno_result", "type": "outfile", 'format': "sequence.profile_table"}  # 注释详细结果表
            ]
        self.add_option(options)
        self.step.add_steps("merge","qs_anno")
        self.result_name = ''

    def check_options(self):
        if not self.option("qs_xml_dir").is_set:
            raise OptionError("Please enter an input file in xml format!", code="31204301")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        self.option('qs_anno_result',os.path.join(self.output_dir,"gene_qs_anno.xls" ))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['gene_qs_anno.xls', 'xls', '序列详细注释文件']
            ])
        super(QsAnnoAgent, self).end()


class QsAnnoTool(Tool):
    def __init__(self, config):
        super(QsAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = self.config.SOFTWARE_DIR + "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + '/annotation/qs_anno.pl'
        self.qs_anno = self.config.SOFTWARE_DIR+ '/database/QS/QS.class'
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.result_name = ''

    def run(self):
        """
        运行
        :return:
        """
        super(QsAnnoTool, self).run()
        self.make_dir()
        self.merge_table()
        self.run_qs_anno()
        self.set_output()
        self.end()

    def merge_table(self):
        self.qs_number = 0
        xml_file = os.listdir(self.option('qs_xml_dir').prop['path'])
        self.result_name = os.path.join(self.work_dir, "qs_align_table.xls")
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        for i in xml_file:
            self.qs_number += 1
            file_path = os.path.join(self.option('qs_xml_dir').prop['path'], i)
            table = xml2table(file_path, self.work_dir + "/tmp_qs_anno/" + "qs_" + str(self.qs_number) + "_table.xls")
            if self.qs_number > 1:
                            os.system("sed -i " +  '1d ' + table)
            cmd = '{} {} {}'.format(self.sh_path, table, self.result_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat" + str(self.qs_number)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i), code="31204301")

    def run_qs_anno(self):
        cmd = '{} {} {} {} {}'.format(self.perl_path, self.perl_script,self.result_name , self.qs_anno,self.work_dir + "/gene_qs_anno.xls")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('running qs_anno end !')
        except subprocess.CalledProcessError:
            self.set_error('running qs_anno error', code="31204302")

    def make_dir(self):
        if os.path.exists(self.work_dir + '/tmp_qs_anno'):
                pass
        else:
                os.mkdir(self.work_dir + '/tmp_qs_anno')

    def set_output(self):
        if os.path.exists(self.output_dir + '/gene_qs_anno.xls'):
            os.remove(self.output_dir + '/gene_qs_anno.xls')
        os.link(self.work_dir + '/gene_qs_anno.xls',self.output_dir + '/gene_qs_anno.xls')
