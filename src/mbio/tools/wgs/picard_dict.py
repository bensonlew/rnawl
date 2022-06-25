# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modify 20180422
# modify 20180514

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class PicardDictAgent(Agent):
    """
    参考基因组的picard建索引
    """
    def __init__(self, parent):
        super(PicardDictAgent, self).__init__(parent)
        options = [
            {"name": "fa", "type": "string"},  # one col
            {"name": "large_genome", "type": "string", "default": "false"}
        ]
        self.add_option(options)
        self.step.add_steps('PicardDict')
    #     self.on('start', self.step_start)
    #     self.on('end', self.step_end)

    # def step_start(self):
    #     self.step.PicardDict.start()
    #     self.step.update()

    # def step_end(self):
    #     self.step.PicardDict.finish()
    #     self.step.update()

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '50G'    # 设置为5 10 20G运行被迫终止！！

    def end(self):
        super(PicardDictAgent, self).end()


class PicardDictTool(Tool):
    def __init__(self, config):
        super(PicardDictTool, self).__init__(config)
        self.java_path = "program/sun_jdk1.8.0/bin/java"
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/picard.jar"
        self.samtools_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/samtools"

    def picarddict(self):
        """
        :return:
        """
        cmd = "{} -jar {} CreateSequenceDictionary REFERENCE={}" \
            .format(self.java_path, self.picard_path, self.option("fa"),
                    self.output_dir)
        self.logger.info(cmd)
        self.logger.info("开始进行PicardDict")
        command = self.add_command("picarddict", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("PicardDict完成！")
        else:
            if self.check_rerun():
                self.logger.info("返回码为1，并且检测是大基因组，使用samtools dict建索引")
                self.samtoolsdict()
            else:
                self.set_error("PicardDict出错！", code="34504201")

    def samtoolsdict(self):
        """
        使用samtools去进行建立dict文件，一般是针对大参考基因组
        modified by hongdong@20200323
        samtools dict ref.fa -o ref.dict
        :return:
        """
        cmd = "{} dict {} -o {}".format(self.samtools_path, self.option("fa"),
                                        os.path.join(os.path.dirname(self.option("fa")), 'ref.dict'))
        self.logger.info(cmd)
        self.logger.info("开始进行PicardDict")
        command = self.add_command("samtoolsdict", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("samtoolsdict完成！")
        else:
            self.set_error("PicardDict出错！", code="34504201")

    def check_rerun(self):
        """
        "Exception in thread "main" java.lang.NegativeArraySizeException"就重新运行一次
        :return:
        """
        with open(self.work_dir + "/picarddict.o", 'r') as r:
            data = r.readlines()
            for line in data:
                if re.match(r".*java\.lang\.NegativeArraySizeException.*", line):
                    self.logger.info("error：{}".format(line))
                    return True
        return False

    def run(self):
        super(PicardDictTool, self).run()
        self.picarddict()
        self.end()
