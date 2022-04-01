# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.06.01

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class CanuAssembleAgent(Agent):
    """
    真菌基因三代数据subreads.fq文件canu软件组装
    """
    def __init__(self, parent):
        super(CanuAssembleAgent, self).__init__(parent)
        options = [
            {"name": "subread_fq", "type": "infile", "format": "sequence.fastq"},  ##三代数据subreads.fq文件
            {"name": "genomeSize", "type": "string"},  ##基因组大小
            {"name": "sample_name", "type": "string"},
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("subread_fq").is_set:
            raise OptionError("必须添加三代数据subreads.fq文件！", code="32100401")
        if not self.option("genomeSize"):
            raise OptionError("必须添加genomeSize的参数，基因组大小！", code="32100402")

    def set_resource(self):
        self._cpu = 4
        memory = 80
        fq = self.option('subread_fq').prop['path']
        fq_size = os.path.getsize(fq)
        if fq_size > 10000000000:
            memory += (fq_size - 10000000000)/1000000000  #超过10G ，每增加1G ，内存加1G
        self._memory = str(memory) + 'G'


    def end(self):
        super(CanuAssembleAgent, self).end()

class CanuAssembleTool(Tool):
    def __init__(self, config):
        super(CanuAssembleTool, self).__init__(config)
        self.subread_fq = self.option("subread_fq").prop['path']
        if re.search("M", self.option('genomeSize')):
            self.genomeSize = self.option('genomeSize')
        else:
            self.genomeSize = self.option("genomeSize")+"M"
        self.sample_name = self.option("sample_name")
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/canu-1.7/Linux-amd64/bin')
        self.set_environ(CLASSPATH='.:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/dt.jar:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/tools.jar')
        self.canu = "bioinfo/Genomic/Sofware/canu-1.7/Linux-amd64/bin/canu"

    def run_canu_correct(self):
        cmd = '{} -correct -p {} -d assemble genomeSize={} -pacbio-raw {} gnuplotTested=true useGrid=false'.format(self.canu,self.sample_name,self.genomeSize,self.subread_fq)
        command = self.add_command("run_canu_correct", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_canu_correct运行完成")
        else:
            self.set_error("run_canu_correct运行出错!", code="32100401")

    def run_canu_trim(self):
        cmd = '{} -trim -p {} -d assemble genomeSize={} -pacbio-corrected assemble/{}.correctedReads.fasta.gz gnuplotTested=true useGrid=false'.format(self.canu,self.sample_name,self.genomeSize,self.sample_name)
        command = self.add_command("run_canu_trim", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_canu_trim运行完成")
        else:
            self.set_error("run_canu_trim运行出错!", code="32100402")

    def run_canu_assemble(self):
        cmd = '{} -assemble -p {} -d assemble-erate-0.013 genomeSize={} correctedErrorRate=0.013 -pacbio-corrected assemble/{}.trimmedReads.fasta.gz gnuplotTested=true useGrid=false'.format(self.canu,self.sample_name,self.genomeSize,self.sample_name)
        command = self.add_command("run_canu_assemble", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_canu_assemble运行完成")
        else:
            self.set_error("run_canu_assemble运行出错!", code="32100403")

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + self.sample_name + '.contigs.fasta'):
            os.remove(self.output_dir + '/' + self.sample_name + '.contigs.fasta')
        os.link(self.work_dir + '/assemble-erate-0.013/' + self.sample_name + '.contigs.fasta',self.output_dir + '/' + self.sample_name + '.contigs.fasta')

    def run(self):
        super(CanuAssembleTool, self).run()
        self.run_canu_correct()
        self.run_canu_trim()
        self.run_canu_assemble()
        self.set_output()
        self.end()