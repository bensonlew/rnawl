# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import subprocess
import re


class MetaSourcetrackerAgent(Agent):
	"""
	meta之微生物组成来源比例分析
	version v1
	author: zhouxuan
	last_modify: 2016.12.19
	"""

	def __init__(self, parent):
		super(MetaSourcetrackerAgent, self).__init__(parent)
		options = [
			# {"name": "otu_table", "type": "infile", "format": "meta.beta_diversity.qiime_table"},  # 输入的OTU文件
			{"name": "otu_table", "type": "string"},  # 输入的OTU文件,在workflow里面处理过的OTU表格路径
			{"name": "map_table", "type": "infile", "format": "meta.otu.group_table"},  # 输入的map文件
			{"name": "s", "type": "string"}  # 物种筛选系数
		]
		self.add_option(options)
		self.step.add_steps("metasourcetracker")
		self.on('start', self.stepstart)
		self.on('end', self.stepfinish)

	def stepstart(self):
		self.step.metasourcetracker.start()
		self.step.update()

	def stepfinish(self):
		self.step.metasourcetracker.finish()
		self.step.update()

	def check_options(self):
		"""
		重写参数检测函数
		:return:
		"""
		if not self.option('otu_table'):
			raise OptionError('必须输入正确的物种（或者otu）信息表', code="32702301")
		if not self.option('map_table'):
			raise OptionError('必须输入map信息表', code="32702302")
		return True

	def set_resource(self):
		"""
		设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
		:return:
		"""
		self._cpu = 10
		self._memory = "20G"

	def end(self):
		result_dir = self.add_upload_dir(self.output_dir)
		result_dir.add_relpath_rules([
			[".", "", "结果输出目录"],
		])
		super(MetaSourcetrackerAgent, self).end()


class MetaSourcetrackerTool(Tool):
	def __init__(self, config):
		super(MetaSourcetrackerTool, self).__init__(config)
		self._version = "v1"
		self.biom_path = 'miniconda2/bin/biom'  #  /mnt/ilustre/users/sanger-dev/app/ app/program/R-3.3.1/lib64/R/bin/exec/R
		self.r_path = 'program/R-3.3.1/bin/Rscript'
		self.filter_otus_from_otu_table_path = 'miniconda2/bin/filter_otus_from_otu_table.py'
		self.sourcetracker_for_qiime_path = self.config.SOFTWARE_DIR + '/bioinfo/meta/sourcetracker-1.0.0-release/sourcetracker_for_qiime.r'

	def run(self):
		"""
		运行
		:return:
		"""
		super(MetaSourcetrackerTool, self).run()
		self.run_metasourcetracker()
		# self.set_output()
		self.end()

	def run_metasourcetracker(self):
		"""
		运行相关脚本和软件，进行微生物组成来源比例分析
		"""
		if self.option('s') != "0" and self.option('s') != "":
			cmd1 = self.biom_path + (' convert -i %s -o temp.biom --table-type "OTU table" --to-hdf5' % (self.option('otu_table')))
			self.logger.info('运行biom，进行OTU转biom')
			command1 = self.add_command("convert1_cmd", cmd1).run()
			self.wait(command1)
			if command1.return_code == 0:
				self.logger.info("OTU转biom运行完成!")
			else:
				self.set_error("OTU转biom运行出错!", code="32702303")
				raise Exception("OTU转biom运行出错!")
			cmd2 = self.filter_otus_from_otu_table_path + (' -i temp.biom -o filtered.biom -s %s' % (self.option('s')))
			self.logger.info('python脚本，进行物种筛选')
			command2 = self.add_command("python_cmd", cmd2).run()
			self.wait(command2)
			if command2.return_code == 0:
				self.logger.info("物种筛选运行完成!")
			else:
				self.set_error("物种筛选运行出错!", code="32702304")
				raise Exception("物种筛选运行出错!")
			cmd3 = self.biom_path + (' convert -i filtered.biom -o filtered.txt --table-type "OTU table" --to-tsv')
			self.logger.info('运行biom，进行biom转OTU')
			command3 = self.add_command("convert2_cmd", cmd3).run()
			self.wait(command3)
			if command3.return_code == 0:
				self.logger.info("biom转OTU运行完成!")
			else:
				self.set_error("biom转OTU运行出错!", code="32702305")
				raise Exception("biom转OTU运行出错!")
			cmd4 = self.r_path + (' %s -i filtered.txt -m %s -o %s' %
								 (self.sourcetracker_for_qiime_path,
								  self.option('map_table').prop['path'], self.output_dir))
		else:
			cmd4 = self.r_path + (
				' %s -i %s -m %s -o %s' % (self.sourcetracker_for_qiime_path,
					self.option('otu_table'), self.option('map_table').prop['path'], self.output_dir))
		self.logger.info('运行R，进行结果数据生成')
		command4 = self.add_command("r_cmd", cmd4).run()
		self.wait(command4)
		if command4.return_code == 0:
			self.logger.info("结果数据生成运行完成!")
		else:
			self.set_error("结果数据生成运行出错!", code="32702306")
			raise Exception("结果数据生成运行出错!")