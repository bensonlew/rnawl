# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.3.8

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class PacbioQcAgent(Agent):
    """
    细菌基因组三代数据质控tool
    """
    def __init__(self, parent):
        super(PacbioQcAgent, self).__init__(parent)
        options = [
            {"name": "input_fofn", "type": "infile", "format": "assembly.fofn"},  # 输入文件,每个样品三代h5文件list
            {"name": "sample_name", "type": "string"}  # 样品名称
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_fofn").is_set:
            raise OptionError("必须添加样品的list文件！")
        if not self.option("sample_name"):
            raise OptionError("必须添加样品名称！")

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(PacbioQcAgent, self).end()

class PacbioQcTool(Tool):
    def __init__(self, config):
        super(PacbioQcTool, self).__init__(config)
        self.fofn = self.option("input_fofn").prop['path']
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/bacgenome/'
        self.set_env = self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/smrtanalysis/smrtanalysis/current/etc/setup.sh'
        os.system("source %s" % self.set_env)
        self.sample_name = self.option("sample_name")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.fofnToSmrt= "/bioinfo/Genomic/Sofware/smrtanalysis/smrtanalysis/current/analysis/bin/fofnToSmrtpipeInput.py"
        self.smrtpipe = "/bioinfo/Genomic/Sofware/smrtanalysis/smrtanalysis/current/analysis/bin/smrtpipe.py"
        self.seqstat = "/bioinfo/Genomic/Sofware/seqstat"
        self.xml_file = self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/smrtanalysis/settings.xml'
        self.otu = self.work_dir + '/input.xml'
        self.out_file = self.work_dir + '/pacbio_qc'
        self.pacbio_summary = self.work_dir + '/pacbio.stat.xls'

    def run_fofn_to_xml(self):
        cmd = self.sh_path + 'fofn_to_xml.sh' + ' ' + self.config.SOFTWARE_DIR +self.fofnToSmrt + ' ' + self.fofn + ' ' + self.otu
        command = self.add_command("run_fofn_to_xml", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_fofn_to_xml运行完成")
        else:
            self.set_error("run_fofn_to_xml运行出错!")

    def run_smrtpipe(self):
        os.system("cp %s settings.xml" % self.xml_file)
        os.system("mkdir %s" % self.out_file)
        cmd = self.sh_path + 'smrtpipe.sh' + ' ' + self.config.SOFTWARE_DIR + self.smrtpipe + ' ' + self.out_file
        #cmd = '{} --distribute -D NPROC=14 -D CLUSTER=BASH -D MAX_THREADS=12 --output={} --params=settings.xml xml:input.xml'.format(self.smrtpipe,self.out_file)
        command = self.add_command("run_smrtpipe", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_smrtpipe运行完成")
        else:
            self.set_error("run_smrtpipe运行出错!")

    def run_summary(self):
        cmd = self.sh_path + 'pacbio_stat.sh' + ' ' + self.config.SOFTWARE_DIR +self.seqstat + ' ' + self.out_file + '/data/filtered_subreads.fasta' + ' ' + self.pacbio_summary
        command = self.add_command("run_summary", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_summary运行完成")
        else:
            self.set_error("run_summary运行出错!")

    def run_summary2(self):
        cmd = '{} {}summary_stat.pl {} {}'.format(self.perl_path,self.perl_script,self.pacbio_summary,self.sample_name + '_PacBio_statistics.xls')
        command = self.add_command("run_summary2", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_summary2运行完成")
        else:
            self.set_error("run_summary2运行出错!")

    def run_len_graphic(self):
        cmd = '{} {}pacbio_length.pl {} {} 100 100 3 4 pacbio'.format(self.perl_path,self.perl_script,self.out_file + '/data/filtered_summary.csv',self.out_file + '/data/filtered_subread_summary.csv')
        command = self.add_command("run_len_graphic", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_len_graphic运行完成")
        else:
            self.set_error("run_len_graphic运行出错!")

    def run_qv_graphic(self):
        cmd = '{} {}pacbio_QV.pl {} 0.001 0.001 4 3 pacbio'.format(self.perl_path,self.perl_script,self.out_file + '/data/filtered_summary.csv')
        command = self.add_command("run_qv_graphic", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_qv_graphic运行完成")
        else:
            self.set_error("run_qv_graphic运行出错!")

    def set_output(self):
        for file in ['pacbio.raw.qv.xls','pacbio.raw.len.xls','pacbio.clean.qv.xls','pacbio.clean.len.xls',self.sample_name + '_PacBio_statistics.xls']:
            if os.path.exists(self.output_dir + '/' + file):
                os.remove(self.output_dir + '/' + file)
            os.link(self.work_dir + '/' + file,self.output_dir + '/' + file)

    def run(self):
        super(PacbioQcTool, self).run()
        self.run_fofn_to_xml()
        self.run_smrtpipe()
        self.run_summary()
        self.run_summary2()
        self.run_len_graphic()
        self.run_qv_graphic()
        self.set_output()
        self.end()