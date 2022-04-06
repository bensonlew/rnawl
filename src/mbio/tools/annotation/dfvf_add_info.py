# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
import shutil
import xml.etree.ElementTree as ET
from biocluster.config import Config
from biocluster.core.exceptions import OptionError


class DfvfAddInfoAgent(Agent):
    """
    author: guanqing.zou
    last_modify: 20180528
    """
    def __init__(self, parent):
        super(DfvfAddInfoAgent, self).__init__(parent)
        options = [
            # {"name": "bsnout", "type": "infile", },
            {"name": "bsnout", "type": "string", "default": ""},
            # {"name": "database", "type": "infile", }
            #{"name": "database", "type": "string", "default": ""}
            ]
        self.add_option(options)
        self.step.add_steps('addinfo')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        #self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def step_start(self):
        self.step.addinfo.start()
        self.step.update()

    def step_end(self):
        self.step.addinfo.finish()
        self.step.update()

    def check_options(self):
        #if not self.option("bsnout").is_set:
        if not self.option("bsnout"):
            raise OptionError("必须设置参数bsnout", code="31201401")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '1G'

    def end(self):
        super(DfvfAddInfoAgent, self).end()


class DfvfAddInfoTool(Tool):
    def __init__(self, config):
        super(DfvfAddInfoTool, self).__init__(config)
        self.data_path = os.path.join(self.config.SOFTWARE_DIR, "database/DFVF/DFVF.data")
        self.python_path = os.path.join(self.config.SOFTWARE_DIR, "/miniconda2/bin/python")

    def run_addinfo(self):
        outputfile = os.path.join(self.output_dir, "result.out.xls")
        cmd_path = os.path.join(self.config.PACKAGE_DIR, "annotation/dfvf_merge_info.py")
        cmd_str = "{} {} {} {} {}".format(self.python_path,cmd_path, self.data_path, self.option("bsnout"), outputfile)

        self.logger.info("开始运行{}".format(cmd_str))
        command = self.add_command("addinfo",cmd_str)
        command.run()
        self.wait(command)
        if command.return_code == 0:
                self.logger.info("{} 运行完成".format(cmd_path))
                self.logger.info(outputfile)
        else:
                self.set_error("%s 运行失败", variables=(cmd_path), code="31201401")


    def run(self):
        """
        运行
        :return:
        """
        super(DfvfAddInfoTool, self).run()
        self.run_addinfo()
        self.end()


