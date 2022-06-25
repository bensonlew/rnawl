# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.06.01

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class CanuAssemAgent(Agent):
    """
    细菌基因三代数据subreads.fq文件canu软件组装
    """
    def __init__(self, parent):
        super(CanuAssemAgent, self).__init__(parent)
        options = [
            {"name": "data_type", "type": "string"},
            {"name": "subread_fq", "type": "infile", "format": "sequence.fastq"},  ##三代数据subreads.fq文件
            {"name": "genomeSize", "type": "string"},  ##基因组大小
            {"name": "sample_name", "type": "string"},
            {"name": "error_rate", "type": "string", "default": "0.025"},  # canu拼接的错误率参数
            {"name": "cor_min_coverage", "type": "string", "default": "default",
             'choose': ['default', '0', '1', '2', "3", "4"]},  # canu param
            {"name": "cor_mhap_sensitivity", "type": "string", "default": "default",
             "choose": ["default", "low", "normal", "high"]}  # canu param
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("subread_fq").is_set:
            raise OptionError("必须添加三代数据subreads.fq文件！", code="32100401")
        if not self.option("genomeSize"):
            raise OptionError("必须添加genomeSize的参数，基因组大小！", code="32100402")

    def set_resource(self):
        self._cpu = 4
        memory = 50
        fq = self.option('subread_fq').prop['path']
        fq_size = os.path.getsize(fq)
        if fq_size > 5000000000:
            memory += (fq_size - 5000000000)/1000000000*10  #超过5G ，每增加1G ，内存加10G
        self._memory = str(memory) + 'G'


    def end(self):
        super(CanuAssemAgent, self).end()

class CanuAssemTool(Tool):
    def __init__(self, config):
        super(CanuAssemTool, self).__init__(config)
        self.subread_fq = self.option("subread_fq").prop['path']
        self.genomeSize = self.option("genomeSize")
        self.sample_name = self.option("sample_name")
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/bin:'+ self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/gnuplot-5.2.3/bin:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/canu-2.0/Linux-amd64/bin:' + self.config.SOFTWARE_DIR + '/miniconda2/bin')
        self.set_environ(CLASSPATH='.:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/dt.jar:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/tools.jar')
        self.set_environ(PERL5LIB=self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/lib/site_perl/5.24.0')
        self.set_environ(LIBRARY_PATH= self.config.SOFTWARE_DIR + '/gcc/7.2.0/lib64')
        self.canu = "bioinfo/Genomic/Sofware/canu-2.0/Linux-amd64/bin/canu"


    def run_canu_correct(self):
        cmd = '{} -correct -p {} -d assemble genomeSize={} useGrid=false'.format(self.canu, self.sample_name, self.genomeSize)
        if self.option("data_type") in ["pacbio", "Pacbio"]:
            cmd += " -pacbio {}".format(self.subread_fq)
        elif self.option("data_type") in ["Nanopore", "nanopore"]:
            cmd += " -nanopore {}".format(self.subread_fq)
        command = self.add_command("run_canu_correct", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_canu_correct运行完成")
        else:
            self.set_error("run_canu_correct运行出错!", code="32100401")

    def run_canu_trim(self):
        cmd = '{} -trim -p {} -d assemble genomeSize={} useGrid=false'.format(self.canu,self.sample_name, self.genomeSize)
        if self.option("data_type") in ["pacbio", "Pacbio"]:
            cmd += " -pacbio-corrected assemble/{}.correctedReads.fasta.gz".format(self.sample_name)
        elif self.option("data_type") in ["Nanopore", "nanopore"]:
            cmd += " -nanopore-corrected assemble/{}.correctedReads.fasta.gz".format(self.sample_name)
        command = self.add_command("run_canu_trim", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_canu_trim运行完成")
        else:
            self.set_error("run_canu_trim运行出错!", code="32100402")


    def run_canu_assemble(self):
        cmd = '{} -assemble -p {} -d assemble-erate-{} genomeSize={} useGrid=false'.format(self.canu,self.sample_name,str(self.option("error_rate")),self.genomeSize, str(self.option("error_rate")))
        if self.option("data_type") in ["pacbio", "Pacbio"]:
            cmd += " -pacbio-corrected assemble/{}.trimmedReads.fasta.gz ".format(self.sample_name)
        elif self.option("data_type") in ["Nanopore", "nanopore"]:
            cmd += " -nanopore-corrected assemble/{}.trimmedReads.fasta.gz ".format(self.sample_name)
        if self.option("cor_min_coverage") != "default":
            cmd += " corMinCoverage={} ".format(self.option("cor_min_coverage"))
        if self.option("cor_mhap_sensitivity") != "default":
            cmd += " corMhapSensitivity={} ".format(self.option("cor_mhap_sensitivity"))
        command = self.add_command("run_canu_assemble", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_canu_assemble运行完成")
        else:
            self.set_error("run_canu_assemble运行出错!", code="32100403")

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + self.sample_name + '.contigs.fasta'):
            os.remove(self.output_dir + '/' + self.sample_name + '.contigs.fasta')
        os.link(self.work_dir + '/assemble-erate-{}/'.format(str(self.option("error_rate"))) + self.sample_name + '.contigs.fasta',self.output_dir + '/' + self.sample_name + '.scaffold.fna')

    def run(self):
        super(CanuAssemTool, self).run()
        self.run_canu_correct()
        self.run_canu_trim()
        self.run_canu_assemble()
        self.set_output()
        self.end()