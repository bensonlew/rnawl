# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.18

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class GihunterPredictAgent(Agent):
    """
    根据输入文件为组装序列,gff文件，rRNA文件预测基因组岛island
    """

    def __init__(self, parent):
        super(GihunterPredictAgent, self).__init__(parent)
        options = [
            {"name": "fna", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "ptt", "type": "string"},  #类似gff
            {"name": "rnt", "type": "string", },   #rNA的统计文件
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fna").is_set:
            raise OptionError("必须设置参数fna文件!")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(GihunterPredictAgent, self).end()


class GihunterPredictTool(Tool):
    def __init__(self, config):
        super(GihunterPredictTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin:' + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/GIHunter/alien_hunter-1.7:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/GIHunter:" + self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "program/Python35/bin"
        self.lib = self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64:"
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0')
        self.set_environ(CLASSPATH='.:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/dt.jar:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/tools.jar')
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.gihunter = "/bioinfo/Genomic/mobile_genetic_elements/GIHunter/GIHunter"
        self.fasta = os.path.basename(self.option("fna").prop['path'])
        self.ptt = os.path.basename(self.option("ptt"))
        self.rnt = os.path.basename(self.option("rnt"))
        self.logger.info(self.fasta)
        self.sample = os.path.basename(self.option("fna").prop['path']).split(".fna")[0]

    def run_gihunter(self):
        if os.path.exists(self.work_dir + "/Genome_inputs"):
            shutil.rmtree(self.work_dir + "/Genome_inputs")
        os.mkdir(self.work_dir + "/Genome_inputs")
        for i in [self.option("fna").prop['path'], self.option("ptt"), self.option("rnt")]:
            self.logger.info(i)
            os.link(i, self.work_dir + "/Genome_inputs/" + os.path.basename(i))
        cmd = "{} {} {} {} {} ".format(self.gihunter, self.work_dir + "/Genome_inputs/" + self.fasta, self.work_dir + "/Genome_inputs/" + self.ptt, self.work_dir + "/Genome_inputs/" + self.rnt, self.sample)
        command = self.add_command("run_gihunter", cmd).run()
        self.logger.info(cmd)
        self.wait(command)
        if command.return_code == 0:
            self.run_summary()
            self.logger.info("run_gihunter运行完成！")
        else:
            self.set_error("run_gihunter运行完成运行出错!")

    def run_summary(self):
        if os.path.exists(self.work_dir + "/Genome_outputs/" + self.sample + "/" + self.sample + "_GIs.txt") >0:
            with open(self.work_dir + "/Genome_outputs/" + self.sample + "/" + self.sample + "_GIs.txt", "r") as f, open(self.work_dir + "/" + self.sample + "_GIs.txt", "w") as g:
                lines = f.readlines()
                for line in lines:
                    lin = line.strip().split("\t")
                    g.write("{}\t{}\t{}\t{}\n".format(self.sample, "Gihunter", lin[0], lin[1]))

    def set_output(self):
        if os.path.exists(self.output_dir + "/" + self.sample + "_GIs.txt"):
            os.remove(self.output_dir + "/" + self.sample + "_GIs.txt")
        if os.path.exists(self.work_dir + "/" + self.sample + "_GIs.txt"):
            os.link(self.work_dir + "/" + self.sample + "_GIs.txt", self.output_dir + "/" + self.sample + "_GIs.txt")
        self.option("out",self.output_dir + "/" + self.sample + "_GIs.txt")

    def run(self):
        super(GihunterPredictTool, self).run()
        self.run_gihunter()
        self.set_output()
        self.end()