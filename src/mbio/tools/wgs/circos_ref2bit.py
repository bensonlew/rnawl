# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180517

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class CircosRef2bitAgent(Agent):
    """
    wgs里基因组浏览器需要这个文件
    """
    def __init__(self,parent):
        super(CircosRef2bitAgent,self).__init__(parent)
        options = [
            {"name": "reffa", "type": "string"},  
        ]
        self.add_option(options)
        self.step.add_steps('CircosRef2bit')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.CircosRef2bit.start()
        self.step.update()

    def step_end(self):
        self.step.CircosRef2bit.finish()
        self.step.update()

    def check_options(self):
        if not self.option("reffa"):
            raise OptionError("请设置reffa", code="34501301") # 必须有

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(CircosRef2bitAgent,self).end()
########

class CircosRef2bitTool(Tool):
    def __init__(self, config):
        super(CircosRef2bitTool, self).__init__(config)
        self.ref2bit_path = '/bioinfo/align/ucsc_tools/faToTwoBit'  # self.config.SOFTWARE_DIR 就是到app/的路径

    def CircosRef2bit(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{} {} {}"\
            .format(self.ref2bit_path, self.option("reffa"), self.output_dir + "/ref.2bit")
        self.logger.info(cmd)
        self.logger.info("开始进行CircosRef2bit")
        command = self.add_command("circosref2bit", cmd).run()  # nranno必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info(cmd)
            self.logger.info("CircosRef2bit完成！")
        else:
            self.set_error("CircosRef2bit出错！", code="34501301")
            self.set_error("CircosRef2bit出错！", code="34501304")

    def run(self):
        super(CircosRef2bitTool, self).run()
        self.CircosRef2bit()
        self.end()
