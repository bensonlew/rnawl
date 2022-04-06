# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/5/6'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
import json


class GetFaFromJsonAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(GetFaFromJsonAgent, self).__init__(parent)
        options = [
            {"name": "chr_fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "genome_json", "type": "string", "required": True},
            {"name": "log", "type": "infile", "format": "sequence.profile_table", "required": True},
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "out_log", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class GetFaFromJsonTool(Tool):
    def __init__(self, config):
        super(GetFaFromJsonTool, self).__init__(config)
        self.script = self.config.PACKAGE_DIR + '/bacgenome/get_fa_by_json.py'
        self.python_path = '/miniconda2/bin/python'
        self.json = json.dumps(self.option("genome_json"))

    def run_fa(self):
        """
        description
        :return:
        """
        if self.option("log").is_set:
            cmd = self.python_path + " " + self.script + " -fa %s -json %s -log %s -seq_prefix Scaffold -o_prefix %s" % (self.option("chr_fa").prop["path"], self.json, self.option("log").prop["path"], self.output_dir + "/chr")
        else:
            cmd = self.python_path + " " + self.script + " -fa %s -json %s -seq_prefix Scaffold -o_prefix %s" % (
            self.option("chr_fa").prop["path"], self.json, self.output_dir + "/chr")
        if self.config.DBVersion:
            cmd += ' -mongo ' + str(self.config.DBVersion)
        command = self.add_command("get_fa", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fa筛选完成")
        else:
            self.set_error("fa筛选出错")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.option("out_fa").set_path(self.output_dir + "/chr.fa")
        self.option("out_log").set_path(self.output_dir + "/chr.log")

    def run(self):
        super(GetFaFromJsonTool, self).run()
        self.run_fa()
        self.set_output()
        self.end()