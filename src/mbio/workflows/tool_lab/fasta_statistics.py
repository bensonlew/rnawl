# -*- coding: utf-8 -*-
# __author__ = 'zzg'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class FastaStatisticsWorkflow(Workflow):
	"""
	����ͳ�ƹ�����

	"""
	def __init__(self, wsheet_object):
		self._sheet = wsheet_object
		super(FastaStatisticsWorkflow, self).__init__(wsheet_object)
		options = [
			{"name": "input_fa", "type": "infile","format": "sequence.fasta"},  # �����ļ�������ֵ
			{"name": "main_id", "type": "string"},
			{'name': 'update_info', 'type': 'string'},
			{"name": "task_id", "type": "string"}
		]
		self.add_option(options)
		self.revise_infiles()
		self.set_options(self._sheet.options())
		self.FastaStatistics = self.add_tool("tool_lab.fasta_statistics")

	def check_options(self):	
		return True
	
	def run(self):
		self.run_FastaStatistics()
		super(FastaStatisticsWorkflow, self).run()
		
	def run_FastaStatistics(self):
		opts = {
			"input_fa": self.option("input_fa")
		}
		self.FastaStatistics.set_options(opts)
		self.FastaStatistics.on("end", self.set_output)
		self.FastaStatistics.run()

	def set_output(self):
		for file in os.listdir(self.FastaStatistics.output_dir):
			os.link(os.path.join(self.FastaStatistics.output_dir, file), os.path.join(self.output_dir, file))
		self.set_db()

	def set_db(self):
		self.logger.info("��ʼ����")
		api_FastaStatistics = self.api.api("tool_lab.fasta_statistics")
		api_FastaStatistics.add_FastaStatistics(self.output_dir + '/out.xlsx',main_id = self.option('main_id'))
		self.logger.info("�������")
		self.end()

	def end(self):
		result_dir = self.add_upload_dir(self.output_dir)
		super(FastaStatisticsWorkflow, self).end()