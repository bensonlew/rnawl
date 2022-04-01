# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.1.2
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class CrisprcasAgent(Agent):
    """
    细菌基因组Crisprcas预测
    """
    def __init__(self, parent):
        super(CrisprcasAgent, self).__init__(parent)
        options = [
            {"name": "seq", "type": "infile","format": "sequence.fasta"},# 参考序列文件
            {'name': 'sample_name', "type": "string"},  # 样本名
            {"name": "analysis", "type": "string", "default": "uncomplete"}  ###流程分析模式complete，uncomplete
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("seq").is_set:
            raise OptionError("必须添加sequence的序列文件！", code="31401101")

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(CrisprcasAgent, self).end()

class CrisprcasTool(Tool):
    def __init__(self, config):
        super(CrisprcasTool, self).__init__(config)
        self.seq= self.option("seq").prop['path']
        self.sample_name = self.option('sample_name')
        self.set_environ(PATH=self.config.SOFTWARE_DIR +'/program/sun_jdk1.8.0/bin')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.minced = "/bioinfo/Genomic/Sofware/minced-master/minced"
        self.crispras =self.work_dir + "/" + self.sample_name + '.crisprs.xls'

    def run_minced(self):
        cmd = '{} -minNR 2 {} {}'.format(self.minced,self.seq,self.crispras)
        self.logger.info(cmd)
        command = self.add_command("run_minced", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_minced运行完成")
        else:
            self.set_error("run_minced运行出错!", code="31401101")

    def run_result_start(self):
        cmd = '{} {}crisprs-start.pl {} {} {}'.format(self.perl_path,self.perl_script,self.option('analysis'),self.crispras,self.sample_name)
        self.logger.info(cmd)
        command = self.add_command("run_result_start", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_result_start运行完成")
        else:
            self.set_error("run_result_start运行出错!", code="31401102")

    def set_output(self):
        if not os.path.exists(self.output_dir + '/CRISPR_Cas'):
            os.mkdir(self.output_dir + '/CRISPR_Cas')
        gli_path = self.output_dir + "/CRISPR_Cas/" + self.sample_name + "_CRISPR_Cas_summary.xls"
        fre_path = self.output_dir + "/CRISPR_Cas/" + self.sample_name + "_CRISPR_Cas_detail.xls"
        for i in [gli_path,fre_path]:
            if os.path.exists(i):
               os.remove(i)
        os.link(self.work_dir + "/" + self.sample_name + ".summary.xls",gli_path)
        os.link(self.work_dir + "/" + self.sample_name + ".detail.xls",fre_path)

    def run(self):
        super(CrisprcasTool, self).run()
        self.run_minced()
        self.run_result_start()
        self.set_output()
        self.end()