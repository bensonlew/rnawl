# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180828

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class NormalizeFaAgent(Agent):
    """
    用dos2unix和picard处理下载的基因组内部
    """
    def __init__(self, parent):
        super(NormalizeFaAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "large_genome", "type": "string", "default": "false"}
        ]
        self.add_option(options)
        self.step.add_steps('NormalizeFa')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.NormalizeFa.start()
        self.step.update()

    def step_end(self):
        self.step.NormalizeFa.finish()
        self.step.update()

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '50G'

    def end(self):
        super(NormalizeFaAgent, self).end()


class NormalizeFaTool(Tool):
    def __init__(self, config):
        super(NormalizeFaTool, self).__init__(config)
        self.unix_path = "bioinfo/dna_evolution/dos2unix"
        self.java_path = "program/sun_jdk1.8.0/bin/java"
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/picard.jar"

    def run_dos2unix(self):
        """
        unix的fa
        """
        cmd = "{} -n {} {}".format(self.unix_path, self.option('ref_fa').prop['path'],
                                   self.work_dir + "/out.fa")
        self.logger.info(cmd)
        self.logger.info("开始进行run_dos2unix")
        command = self.add_command("dos2unix", cmd).run()   # 必须小写
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_dos2unix完成！")
        else:
            self.set_error("run_dos2unix出错！")

    def run_normalizefa(self):
        """
        picard处理基因组，均一化基因组
        :return:
        """
        cmd = "{} -jar {} NormalizeFasta I={}  O={}" \
            .format(self.java_path, self.picard_path, self.work_dir + "/out.fa",
                    self.output_dir + "/ref.fa")
        self.logger.info(cmd)
        self.logger.info("开始进行NormalizeFa")
        command = self.add_command("normalizefa", cmd, ignore_error=True).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("NormalizeFa完成！")
        else:
            if self.check_rerun():
                self.logger.info("返回码为1，并且检测是大基因组，不再进行NormalizeFasta运算")
                if os.path.exists(self.output_dir + "/ref.fa"):
                    os.remove(self.output_dir + "/ref.fa")
                os.link(self.work_dir + '/out.fa', self.output_dir + "/ref.fa")
            else:
                self.set_error("NormalizeFa出错！")

    def check_rerun(self):
        """
        "Exception in thread "main" java.lang.NegativeArraySizeException"
        :return:
        """
        with open(self.work_dir + "/normalizefa.o", 'r') as r:
            data = r.readlines()
            for line in data:
                if re.match(r".*java\.lang\.NegativeArraySizeException.*", line):
                    self.logger.info("error：{}".format(line))
                    return True
        return False

    def run(self):
        super(NormalizeFaTool, self).run()
        self.run_dos2unix()
        self.run_normalizefa()
        self.end()
