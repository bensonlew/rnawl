# -*_ coding: utf-8 -*-
# __author__ = 'XueQinwen'
# last_modified: 20210907
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os

class PbindexAgent(Agent):
    """
    SMRT Link v10.1
    创建bam.pbi文件
    """
    def __init__(self, parent):
        super(PbindexAgent,self).__init__(parent)
        options = [
            {"name":"input_bam","type":"infile","format":"align.bwa.bam"},
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen
    
    def check_option(self):
        """
        参数检查
        """
        if not self.option("input_bam"):
            raise OptionError("input_bam参数没有填入")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 3
        self._memory = "20G"

    def end(self):
        super(PbindexAgent,self).end()

class PbindexTool(Tool):
    def __init__(self, config):
        super(PbindexTool,self).__init__(config)
        self.pbindex = "program/SmrtLink/smrtlink/smrtcmds/bin/pbindex"

    def run(self):
        super(PbindexTool,self).run()
        self.run_pbindex()
        self.end()

    def run_pbindex(self):
        """
        pbindex input
        """
        self.bam = self.option("input_bam").prop['path']
        cmd = self.pbindex + " " + self.bam
        command = self.add_command("pbindex", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("建索引成功")
        else:
            self.set_error("建索引失败")




