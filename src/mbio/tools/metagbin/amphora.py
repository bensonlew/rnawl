# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.01.14

import os
import re,shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class AmphoraAgent(Agent):
    """
    amphora进行bin的物种注释
    """
    def __init__(self, parent):
        super(AmphoraAgent, self).__init__(parent)
        options = [
            {"name": "bin_fa", "type": "infile", "format": "sequence.fasta"},  # binning的scaffold文件
            {"name": "kingdom", "type": "string", "default": "Bacteria"},  # genome type Kingdom: Bacteria or Archaea
            {"name": "bin_name", "type": "string"},  # bin name
            {"name": "out","type": "outfile","format": "metagbin.file_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bin_fa").is_set:
            raise OptionError("必须设置参数binning的scaffold文件!")
        if self.option("kingdom") == "":
            raise OptionError("必须设置参数kingdom的type：arc or bac!")

    def set_resource(self):
        self._cpu =4
        self._memory = '30G'

    def end(self):
        super(AmphoraAgent, self).end()

class AmphoraTool(Tool):
    def __init__(self, config):
        super(AmphoraTool, self).__init__(config)
        self.path =self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/hmmer-3.0/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/RAxML-master:" + self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/bin:"
        self.set_environ(PATH=self.path,AMPHORA2_home=self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/AMPHORA2-master")
        self.fasta = self.option("bin_fa").prop['path']
        self.perl_path = "/miniconda2/bin/perl"
        self.phyloty = "../../../../../.." + self.config.PACKAGE_DIR + '/metagbin/amphora.sh'
        self.perl_script = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/AMPHORA2-master/Scripts/"
        self.result = self.work_dir + '/' + self.option("bin_name") + '.amphora_anno.xls'
        self.summary =self.config.PACKAGE_DIR + "/metagbin/amphora_anno.pl"
        self.fa_name =os.path.basename(self.fasta)
        if os.path.exists(self.work_dir + '/' + self.fa_name):
            os.remove(self.work_dir + '/' + self.fa_name)
        os.link(self.fasta,self.work_dir + '/' + self.fa_name)

    def run_mrakerscan(self):
        cmd = "{} {}MarkerScanner.pl -DNA ".format(self.perl_path, self.perl_script)
        if self.option("kingdom") in ['Bacteria']:
            cmd += '-Bacteria'
        elif self.option("kingdom") in ['Archaea']:
            cmd += '-Archaea'
        cmd +=" {}".format(self.work_dir + '/' + self.fa_name)
        command = self.add_command("run_mrakerscan", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("amphora比对mraker运行完成！")
        else:
            self.set_error("amphora比对mraker运行完成运行出错!")

    def run_marker_align(self):
        cmd = "{} {}MarkerAlignTrim.pl -WithReference -OutputFormat phylip".format(self.perl_path, self.perl_script)
        command = self.add_command("run_marker_align", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("amphora比对筛选结果运行完成！")
        else:
            self.set_error("amphora比对筛选结果运行完成运行出错!")

    def run_phylotyping(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        os.mkdir(self.work_dir + "/temp")
        cmd = '{} {}Phylotyping.pl -CPUs 4 -tmp_dir {}'.format(self.perl_path, self.perl_script,self.work_dir + "/temp")
        command = self.add_command("run_phylotyping", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            os.system("cat {}*.phylotype >result.phylotype" .format(self.work_dir + "/"))
            self.logger.info("提取bin注释结果文件运行完成!")
        else:
            self.set_error("提取bin注释结果文件运行出错!")

    def run_sum_anno(self):
        cmd = '{} {} {}'.format(self.perl_path, self.summary,self.work_dir)
        command = self.add_command("run_sum_anno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("汇总bin注释结果文件运行完成!")
        else:
            self.set_error("汇总bin注释结果文件运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + self.option("bin_name") + '.amphora_anno.xls'):
            os.remove(self.output_dir + '/' + self.option("bin_name") + '.amphora_anno.xls')
        os.link(self.work_dir + '/bin.anno.xls',self.output_dir + '/' + self.option("bin_name") + '.amphora_anno.xls')
        self.option('out', self.output_dir + '/' + self.option("bin_name") + '.amphora_anno.xls')
    
    def run_re(self):
        if os.path.exists(self.output_dir + '/' + self.option("bin_name") + '.amphora_anno.xls'):
            os.remove(self.output_dir + '/' + self.option("bin_name") + '.amphora_anno.xls')
        os.link(self.work_dir + '/result.phylotype',self.output_dir + '/' + self.option("bin_name") + '.amphora_anno.xls')
        self.option('out', self.output_dir + '/' + self.option("bin_name") + '.amphora_anno.xls')       
          
    def run(self):
        super(AmphoraTool, self).run()
        self.run_mrakerscan()
        self.run_marker_align()
        self.run_phylotyping()
        with open (self.work_dir + "/result.phylotype","r") as f:
            lines = f.readlines()
            if len(lines) <=1:
               self.run_re()
            else:
               self.run_sum_anno()
               self.set_output()
        self.end()
