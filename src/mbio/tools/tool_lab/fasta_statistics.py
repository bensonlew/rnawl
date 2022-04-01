# -*- coding: utf-8 -*-
# __author__ = zhaozhigang
# last_modify:2020.8.4

import re, os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class FastaStatisticsAgent(Agent):
	"""
	统计序列最小长度、最大长度、平均长度、N50、N90、gap率（scaffold）
	"""
	def __init__(self, parent):
		super(FastaStatisticsAgent, self).__init__(parent)
		options = [
			{"name": "input_fa", "type": "infile", "format": "sequence.fasta"},  # 输入文件,序列文件
			{"name": "out_table", "type": "outfile", "format": "sequence.profile_table"}
		]
		self.add_option(options)
		
	def check_options(self):
		pass
	
	def set_resource(self):
		self._cpu = 2
		self._memory = '2G'
"""
	def end(self):
		super(FastaStatisticsAgent, self).end()	
"""	
		
class FastaStatisticsTool(Tool):
	def __init__(self, config):
		super(FastaStatisticsTool, self).__init__(config)
		self.python_path = "program/Python/bin/python"
		self.python_script = self.config.PACKAGE_DIR + "/tool_lab/fasta_information-v3.py"
				
		
	def run_fasta(self):
		cmd = '{} {} -i {} -o {}'.format(self.python_path,self.python_script,self.option("input_fa").prop["path"],self.output_dir + "/out.xlsx")
		self.logger.info(cmd)
		command = self.add_command("run_fasta", cmd).run()
		self.wait(command)
		if command.return_code == 0:
			self.logger.info("start run_fasta")
		else:
			self.set_error("run error")
		
	def set_output(self):
		self.logger.info("set dir")
		self.option("out_table").set_path(self.output_dir + "/out.xlsx")
		self.logger.info("set dir over")
	
	def run(self):
		super(FastaStatisticsTool, self).run()
		self.run_fasta()
		self.set_output()
		self.end()