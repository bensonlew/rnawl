# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.18

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class MgGihunterPredictAgent(Agent):
    """
    根据输入文件为组装序列,gff文件，rRNA文件预测基因组岛island
    """

    def __init__(self, parent):
        super(MgGihunterPredictAgent, self).__init__(parent)
        options = [
            {"name": "dir", "type": "string"},
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("dir"):
            raise OptionError("必须设置参数dir目录!")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(MgGihunterPredictAgent, self).end()


class MgGihunterPredictTool(Tool):
    def __init__(self, config):
        super(MgGihunterPredictTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin:' + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/GIHunter/alien_hunter-1.7:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/GIHunter:" + self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "program/Python35/bin"
        self.lib = self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64:"
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0')
        self.set_environ(CLASSPATH='.:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/dt.jar:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/tools.jar')
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.gihunter = "/bioinfo/Genomic/mobile_genetic_elements/GIHunter/GIHunter"

    def run_gihunter(self):
        if os.path.exists(self.work_dir + "/result"):
            shutil.rmtree(self.work_dir + "/result")
        os.mkdir(self.work_dir + "/result")
        n = 1
        for i in os.listdir(self.option("dir")):
            if os.path.exists(self.work_dir + "/Genome_inputs"):
                shutil.rmtree(self.work_dir + "/Genome_inputs")
            os.mkdir(self.work_dir + "/Genome_inputs")
            for file in os.listdir(self.option("dir")+ "/" + i):
                self.logger.info(file)
                os.link(self.option("dir")+ "/" + i + "/" + file, self.work_dir + "/Genome_inputs/" + file)
            cmd = "{} {} {} {} {} ".format(self.gihunter, self.work_dir + "/Genome_inputs/" + i + '.fna',
                                           self.work_dir + "/Genome_inputs/" + i + ".ptt",
                                           self.work_dir + "/Genome_inputs/" + i + ".rnt", i)
            command = self.add_command("run_gihunter"+str(n), cmd).run()
            self.logger.info(cmd)
            self.wait(command)
            if command.return_code == 0:
                if os.path.getsize(self.work_dir + "/Genome_outputs/" + i + "/" + i + "_GIs.txt") > 0:
                    os.link(self.work_dir + "/Genome_outputs/" + i + "/" + i + "_GIs.txt", self.work_dir + "/result/" + i + "_GIs.txt")
                self.logger.info("run_gihunter运行完成！")
            else:
                self.set_error("run_gihunter运行完成运行出错!")
            shutil.rmtree(self.work_dir + "/Genome_outputs/")
            n += 1


    def run_summary(self):
        with open(self.output_dir + "/gihunter.xls", "w") as g:
            for file in os.listdir(self.work_dir + "/result"):
                location = file.split("_GIs")[0]
                with open(self.work_dir + "/result/" +file, "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        lin = line.strip().split("\t")
                        g.write("{}\t{}\t{}\t{}\n".format(location, "Gihunter", lin[0], lin[1]))

    def set_output(self):
        if os.path.getsize(self.output_dir + "/gihunter.xls") > 0:
            self.option("out",self.output_dir + "/gihunter.xls")

    def run(self):
        super(MgGihunterPredictTool, self).run()
        self.run_gihunter()
        self.run_summary()
        self.set_output()
        self.end()