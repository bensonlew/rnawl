# -*- coding: utf-8 -*-
# __author__ :gaohao

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os,re
import math


class SpringFileAgent(Agent):
    def __init__(self, parent):
        super(SpringFileAgent, self).__init__(parent)
        options = [
            {"name": "input", "type": "string"},
            {"name": "type", "type": "string"},# compress or decompress
        ]
        self.add_option(options)
        self.queue = 'SPRING'

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        return True

    def set_resource(self):  # 后续需要测试确认
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 4
        memory = os.path.getsize(self.option('input'))
        n = memory / (1024 * 1024 * 1024)
        if int(n) <=1:
            self._memory = "5G"
        else:
            if self.option("type") in ["compress"]:
                self._memory = '{}G'.format(int(n)*3)
            elif self.option("type") in ["decompress"]:
                self._memory = '{}G'.format(int(n)*10)



    def end(self):
        super(SpringFileAgent, self).end()


class SpringFileTool(Tool):
    def __init__(self, config):
        super(SpringFileTool, self).__init__(config)
        self.path =self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin"
        self.set_environ(PATH=self.path,LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64")
        self._version = "v1.0"
        self.spring = "/bioinfo/seq/Spring-1.0/build/spring"

    def run(self):
        super(SpringFileTool, self).run()
        if self.option("type") in ["compress"]:
            self.run_compress()
        elif self.option("type") in ["decompress"]:
            self.run_decompress()
        self.end()

    def run_compress(self):
        self.logger.info("运行spring 对文件进行压缩！")
        name = os.path.basename(self.option("input"))
        cmd = "{} -c -i {} -o {} -t 4".format(self.spring, self.option("input"), self.output_dir + "/" + name + ".spr")
        command = self.add_command("run_compress", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_compress运行完成")
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("run_compress运行出错！")

    def run_decompress(self):
        self.logger.info("运行spring 对文件进行解压！")
        name = os.path.basename(self.option("input")).split(".spr")[0]
        cmd = "{} -d -i {} -o {} -t 4".format(self.spring, self.option("input"), self.output_dir + "/" + name)
        command = self.add_command("run_decompress", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_decompress运行完成")
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("run_decompress运行出错！")