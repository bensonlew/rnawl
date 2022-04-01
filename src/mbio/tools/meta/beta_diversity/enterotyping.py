# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import subprocess
import re


class EnterotypingAgent(Agent):
	"""
	meta之样本菌群分型分析
	version v1
	author: zhouxuan
	last_modify: 2016.12.9
	"""

	def __init__(self, parent):
		super(EnterotypingAgent, self).__init__(parent)
		options = [
			{"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的文件
			#{"name": "result_dir", "type": "outfile", "format": "meta.beta_diversity.result_dir"}  # 输出文件夹
		]
		self.add_option(options)
		self.step.add_steps("enterotyping")
		self.on('start', self.stepstart)
		self.on('end', self.stepfinish)

	def stepstart(self):
		self.step.enterotyping.start()
		self.step.update()

	def stepfinish(self):
		self.step.enterotyping.finish()
		self.step.update()

	def check_options(self):
		"""
		重写参数检测函数
		:return:
		"""
		if not self.option('otu_table'):
			raise OptionError('必须输入正确的物种（或者otu）信息表')
		return True

	def set_resource(self):
		"""
		设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
		:return:
		"""
		self._cpu = 1
		self._memory = "5G"

	def end(self):
		result_dir = self.add_upload_dir(self.output_dir)
		result_dir.add_relpath_rules([
			[".", "", "结果输出目录"],
		])
		super(EnterotypingAgent, self).end()


class EnterotypingTool(Tool):
	def __init__(self, config):
		super(EnterotypingTool, self).__init__(config)
		self._version = "v1"
		self.enterotyping_path = '/mnt/ilustre/users/sanger-dev/app/bioinfo/meta/scripts'
		self.perl_path = 'program/perl-5.24.0/bin/perl '
		self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'  # 2016.12.26 by zhouxuan 3
		self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
		self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
		# self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')  # 2016.12.26 by zhouxuan

	def run(self):
		"""
		运行
		:return:
		"""
		super(EnterotypingTool, self).run()
		self.run_enterotyping()
		# self.set_output()
		self.end()

	def run_enterotyping(self):
		"""
		运行perl脚本，进行样本菌群分型分析
		"""
		cmd = self.perl_path + self.enterotyping_path + (
			'/Enterotyping.pl -i %s -o %s' % (self.option('otu_table').prop['path'], self.output_dir))
		self.logger.info('运行perl脚本，进行样本菌群分型分析')
		command = self.add_command("enterotyping_cmd", cmd).run()
		self.wait(command)
		if command.return_code == 0:
			self.logger.info("enterotyping运行完成!")
		else:
			self.set_error("enterotyping运行出错!")

	# def set_output(self):
	# 	"""
	# 	将结果文件复制到output文件夹下面
	# 	:return:
	# 	"""
	# 	self.logger.info("设置结果目录")
	# 	self.option('result_dir').set_path(self.output_dir)
	# 	self.logger.info("设置样本菌群分型分析成功")




