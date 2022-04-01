# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.12

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.mobile_genetic_elements.commom import get_num,link_file


class IntegronPredictAgent(Agent):
    """
    根据输入文件的不同进行分析，进行integron的预测
    """

    def __init__(self, parent):
        super(IntegronPredictAgent, self).__init__(parent)
        options = [
            {"name": "input_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "file_dir", "type": "string"},  # 该文件夹下有*.gff,*.faa,*.fna(基因的gff文件，基因蛋白文件，基因的核酸文件)
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_fa").is_set:
            raise OptionError("必须设置参数input_fa文件!")

    def set_resource(self):
        self._cpu = 4
        self._memory = '10G'

    def end(self):
        super(IntegronPredictAgent, self).end()


class IntegronPredictTool(Tool):
    def __init__(self, config):
        super(IntegronPredictTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/hmmer-3.1b2/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/infernal-1.1.3-linux-intel-gcc/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/Integron_Finder-master/data:"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib:"
        self.set_environ(PATH=self.path,LD_LIBRARY_PATH=self.lib)
        self.python = "/program/Python35/bin/integron_finder"
        self.fasta = self.option("input_fa").prop['path']
        self.names = re.split("\.fna|\.fasta|\.fa",os.path.basename(self.fasta))[0]

    def run_integron(self):
        if os.path.exists(self.work_dir + "/Results_Integron_Finder_" + self.names):
            shutil.rmtree(self.work_dir + "/Results_Integron_Finder_" + self.names)
        if self.option("file_dir"):
            shutil.copytree(self.option("file_dir"), self.work_dir + "/Results_Integron_Finder_" + self.names)
        cmd = "{} --local-max --func-annot --keep-tmp --cpu 4 --outdir {} {}".format(self.python, self.work_dir, self.fasta)
        command = self.add_command("run_integron", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_integron运行完成！")
        else:
            self.set_error("run_integron运行完成运行出错!")

    def set_output(self):
        self.logger.info(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names +".integrons")
        num =get_num(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names +".integrons")
        if num >2:
            link_file(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names + ".integrons", self.output_dir + "/" + self.names + ".integrons")
            link_file(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names + ".summary", self.output_dir + "/" + self.names + ".summary")


    def run(self):
        super(IntegronPredictTool, self).run()
        self.run_integron()
        self.set_output()
        self.end()
