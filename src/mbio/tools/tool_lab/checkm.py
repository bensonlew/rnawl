#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'zzg'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil

class CheckmAgent(Agent):
	"""
	用于评估fasta完整性
	version 1.0
	author: zzg
	last_modify: 2020.9.1
	"""

	def __init__(self, parent):
		super(CheckmAgent, self).__init__(parent)
		options = [
			{"name": "input_file", "type": "infile", "format": "sequence.fasta"},
			{"name": "input_format", "type": "string", "default": "assembly"},
			{"name": "method", "type": "string", "default": "marker set"}
		]
		self.add_option(options)
		self.list =[]


	def check_options(self):
		"""
		检测参数是否正确
		"""
		'''
		if not self.option('genome_dir').is_set and not self.option('protein_dir'):
			raise OptionError("fasta文件夹不存在！")
		if self.option('genome_dir').is_set:
			files = os.listdir(self.option("genome_dir").prop['path'])
			for file in files:
				if file.split(".")[-1] != "fasta":
					raise OptionError("fasta文件后缀必须为.fasta", code="31200102")
		if self.option('protein_dir').is_set:
			files = os.listdir(self.option("genome_dir").prop['path'])
			for file in files:
				if file.split(".")[-1] != "fasta":
					raise OptionError("fasta文件后缀必须为.fasta", code="31200103")
		'''
		return True
		
	def set_resource(self):
		"""
		所需资源
		"""
		self._cpu = 8
		self._memory = '80G'
		
	def end(self):
		super(CheckmAgent, self).end()


class CheckmTool(Tool):
	"""
	version 1.0
	"""
	def __init__(self, config):
		super(CheckmTool, self).__init__(config)
		self.path =self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/CheckM-master:" + self.config.SOFTWARE_DIR + "/program/Python/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/islandpath_dimob/hmmer-3.1b1/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/pplacer-Linux-v1.1.alpha19"
		self.set_environ(PATH=self.path)
		self.checkm = "/program/Python/bin/"
		self.perl = "/program/perl/perls/perl-5.24.0/bin/perl"
		self.perl_script = self.config.PACKAGE_DIR + "/metagbin/checkm_stat.pl"
		self.summary =self.work_dir + '/checkm_summary.xls'
		self.result = self.work_dir + '/result/'
		self.input_dir = os.path.abspath(os.path.join(self.option('input_file').prop['path'], ".."))
		if self.option('method')  in ["marker gene"]:
			self.individual_markers = "--individual_markers"
		else:
			self.individual_markers = ""
			
			

	def run_checkm_genome(self):
		cmd = "{}checkm lineage_wf -x fasta -t 16 --pplacer_threads 16 -f {} {} --tab_table --tmpdir {} {} {}".format(self.checkm,self.summary,self.individual_markers,"./",self.input_dir,self.result)
		self.logger.info(cmd)
		self.logger.info("开始运行run_checkm")
		command = self.add_command("run_checkm_genome", cmd)
		command.run()
		self.wait(command)
		if command.return_code == 0:
			self.logger.info("运行run_checkm完成")
		else:
			self.set_error("运行run_checkm运行出错!")

	def run_checkm_protein(self):
		cmd = "{}checkm lineage_wf -x fasta -t 16 --pplacer_threads 16 -f {} --gene --tab_table {} --tmpdir {} {} {}".format(self.checkm,self.summary,self.individual_markers,'./',self.input_dir,self.result)
		self.logger.info(cmd)
		self.logger.info("开始运行run_checkm")
		command = self.add_command("run_checkm", cmd)
		command.run()
		self.wait(command)
		if command.return_code == 0:
			self.logger.info("运行run_checkm完成")
		else:
			self.set_error("运行run_checkm运行出错!")

	def set_output(self):
		if os.path.exists(self.output_dir + "/" + 'checkm_summary.xls'):
			os.remove(self.output_dir + "/" + 'checkm_summary.xls')
		os.link(self.work_dir + "/" + 'checkm_summary.xls',self.output_dir + "/" + 'checkm_summary.xls')

	
	def run(self):
		"""
		运行
		"""
		super(CheckmTool, self).run()
		if self.option("input_format") in ["protein"]:
			self.run_checkm_protein()
			self.set_output()
			self.end()
		else:
			self.run_checkm_genome()
			self.set_output()
			self.end()