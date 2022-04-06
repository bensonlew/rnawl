# -*- coding: utf-8 -*-
# __author__ = zhaozhigang
# last_modify:2020.8.4

import re, os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class NormalTestAgent(Agent):
	"""
	
	"""
	def __init__(self, parent):
		super(NormalTestAgent, self).__init__(parent)
		options = [
			{"name": "input_data", "type": "infile","format": "tool_lab.simple"},  # 输入文件,样本值
			{"name": "out_table", "type": "outfile","format": "sequence.profile_table"}
		]
		self.add_option(options)
		
	def check_options(self):
		pass
	
	def set_resource(self):
		self._cpu = 2
		self._memory = '10G'
	"""
	def end(self):
		super(NormalTestAgent, self).end()	
	"""	
		
class NormalTestTool(Tool):
	def __init__(self, config):
		super(NormalTestTool, self).__init__(config)
		self.python_path = "miniconda2/bin/python"
		self.python_script = self.config.PACKAGE_DIR + "/tool_lab/ks_testv2.py"

	def run_calccalc(self):
		cmd = '{} {} -i {} -o {}'.format(self.python_path,self.python_script,self.option("input_data").prop["path"],self.output_dir + "/out1.xlsx")
		self.logger.info(cmd)
		command = self.add_command("run_calccalc", cmd).run()
		self.wait(command)
		if command.return_code == 0:
			self.logger.info("start run_calccalc")
		else:
			self.set_error("run error")
		
	def set_output(self):
		self.logger.info("set dir")
		os.link(self.work_dir+"/histogram_X_axis.txt", self.output_dir+"/histogram_X_axis.txt")
		os.link(self.work_dir+"/histogram_Y_axis.txt", self.output_dir+"/histogram_Y_axis.txt")
		os.link(self.work_dir+"/histogram_line_X.txt", self.output_dir+"/histogram_line_X.txt")
		os.link(self.work_dir + "/histogram_line_Y.txt", self.output_dir + "/histogram_line_Y.txt")
		os.link(self.work_dir+"/QQ_line.txt", self.output_dir+"/QQ_line.txt")
		os.link(self.work_dir+"/QQ_X_axis.txt", self.output_dir+"/QQ_X_axis.txt")
		os.link(self.work_dir+"/QQ_Y_axis.txt", self.output_dir+"/QQ_Y_axis.txt")
		self.logger.info("set dir over")
	
	def run(self):
		super(NormalTestTool, self).run()
		self.run_calccalc()
		self.set_output()
		self.end()