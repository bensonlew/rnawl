# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# version 1.0
# last_modify: 2019.04.23

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
            {"name": "read_type", "type": "string", "default": "pacbio-raw"},  ## pacbio-raw nanopore-raw
            {"name": "correctedErrorRate", "type": "string", "default": "0.013"},
            {"name": "genomeSize", "type": "string"},  ##基因组大小
            {"name": "sample_name", "type": "string"},
            {"name": "canu_version", "type": "string", "default": "1.3", 'choose': ["1.3", "1.7"]},  # canu version
            {"name": "corMinCoverage", "type": "string", "default": "default", 'choose': ['default', '0', '1', '2', "3", "4"]},
            {"name": "corMhapSensitivity", "type": "string", "default": "default", "choose": ["default", "low", "normal", "high"]}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("subread_fq").is_set:
            raise OptionError("必须添加三代数据subreads.fq文件！")
        if not self.option("genomeSize"):
            raise OptionError("必须添加genomeSize的参数，基因组大小！")

    def set_resource(self):
        self._cpu = 4
        self._memory = '80G'

    def end(self):
        super(CanuAssembleAgent, self).end()

class CanuAssembleTool(Tool):
    def __init__(self, config):
        super(CanuAssembleTool, self).__init__(config)
        self.logger.debug(os.environ)
        self.subread_fq = self.option("subread_fq").prop['path']
        self.genomeSize = self.option("genomeSize")
        self.error_rate = self.option("correctedErrorRate")
        self.read_type = self.option("read_type")
        self.sample_name = self.option("sample_name")
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/miniconda2/bin")
        self.logger.debug("after set per environ:")  # 计算节点和登录节点的环境变量可能不一致，导致canu依赖的perl和java环境出现问题，此处作为检查
        self.logger.info(os.environ)
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/canu-1.7/Linux-amd64/bin')
        self.set_environ(CLASSPATH='.:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/dt.jar:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/tools.jar')
        self.canu_17 = "bioinfo/Genomic/Sofware/canu-1.7/Linux-amd64/bin/canu"
        self.canu_13 = "bioinfo/Genomic/Sofware/canu-1.3/Linux-amd64/bin/canu"


    def run_canu_assemble(self):
        # cmd = '{}  -p {} -d assemble genomeSize={}m correctedErrorRate={} -{} {} gnuplotTested=true useGrid=false'.format(self.canu_17,self.sample_name,self.genomeSize,self.error_rate, self.read_type, self.option("subread_fq").prop["path"])
        cmd = "{}  -p {} -d assemble genomeSize={}m correctedErrorRate={} gnuplotTested=true useGrid=false".format(self.canu_17,self.sample_name,self.genomeSize,self.error_rate)
        if self.option("corMinCoverage") != "default":
            cmd += " corMinCoverage={} ".format(self.option("corMinCoverage"))
        if self.option("corMhapSensitivity") != "default":
            cmd += " corMhapSensitivity={} ".format(self.option("corMhapSensitivity"))
        cmd += " -{} {}".format(self.read_type, self.option("subread_fq").prop["path"])
        command = self.add_command("run_canu_assemble", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_canu_assemble1.7运行完成")
        else:
            self.set_error("run_canu_assemble1.7运行出错!")

    def run_canu_1_3(self):
        # cmd =  '{}  -p {} -d assemble genomeSize={}m corMinCoverage=2 corMhapSensitivity=high errorRate={} -{} {} useGrid=false'.format(self.canu_13,self.sample_name,self.genomeSize,self.error_rate, self.read_type, self.option("subread_fq").prop["path"])
        cmd = '{} -p {} -d assemble genomeSize={}m errorRate={} useGrid=false '.format(self.canu_13,self.sample_name,self.genomeSize,self.error_rate)
        if self.option("corMinCoverage") != "default":
            cmd += " corMinCoverage={} ".format(self.option("corMinCoverage"))
        if self.option("corMhapSensitivity") != "default":
            cmd += " corMhapSensitivity={} ".format(self.option("corMhapSensitivity"))
        cmd += " -{} {}".format(self.read_type, self.option("subread_fq").prop["path"])
        command = self.add_command("run_canu_assemble", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_canu1.3_assemble运行完成")
        else:
            self.set_error("run_canu1.3_assemble运行出错!")

    def set_output(self):
        if os.path.isfile(self.output_dir + '/' + self.sample_name + '_scaffold.fna'):
            os.remove(self.output_dir + '/' + self.sample_name + '_scaffold.fna')
        os.link(self.work_dir + '/assemble/' + self.sample_name + '.contigs.fasta',self.output_dir + '/' + self.sample_name + '.scaffold.fna')

    def run(self):
        super(CanuAssembleTool, self).run()
        if self.option("canu_version") == "1.3":
            self.run_canu_1_3()
        else:
            self.run_canu_assemble()
        self.set_output()
        self.end()