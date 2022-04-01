# -*_ coding: utf-8 -*-
# __author__ = 'XueQinwen'
# last_modified: 20210907
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os

class CcsAgent(Agent):
    """
    SMRT Link v10.1
    从subreads里提取ccs数据
    """
    def __init__(self, parent):
        super(CcsAgent,self).__init__(parent)
        options = [
            {"name":"input_bam","type":"infile","format":"align.bwa.bam"},
            {"name":"chunk","type":"string"},
            {"name":"ccs_bam","type":"string"}
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen
    
    def check_option(self):
        """
        参数检查
        """
        if not self.option("input_bam"):
            raise OptionError("input_bam参数没有填入")
        if not self.option("chunk"):
            raise OptionError("chunk参数没有填入")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 15
        self._memory = "40G"

    def end(self):
        super(CcsAgent,self).end()

class CcsTool(Tool):
    def __init__(self, config):
        super(CcsTool,self).__init__(config)
        self.ccs = "program/SmrtLink/smrtlink/smrtcmds/bin/ccs"

    def run(self):
        super(CcsTool,self).run()
        self.run_ccs()
        self.end()

    def run_ccs(self):
        """
        ccs -j 15 --min-length 100 --max-length 100 --max-length 100000 --min-passes 3 --min-rq 0.99 
        subreads.bam  ccs.bam --report-file ccs.txt --chunk n/20
        """
        self.ccs_bam = self.option("ccs_bam")
        cmd = '{} -j 15 --min-length 100 --max-length 100 --max-length 100000 --min-passes 3 --min-rq 0.99 '.format(self.ccs)
        cmd += ' {} {} --report-file {}/ccs.txt --chunk {}/10'.format(self.option("input_bam").prop['path'],self.ccs_bam,self.work_dir,self.option('chunk'))
        command = self.add_command("ccs{}".format(self.option('chunk')),cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("从subread中提取ccs完成")
        else:
            self.set_error("提取ccs数据发生问题")




