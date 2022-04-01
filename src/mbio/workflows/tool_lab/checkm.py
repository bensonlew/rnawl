# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import shutil
from bson.objectid import ObjectId


class CheckmWorkflow(Workflow):
	"""
	checkM

	"""
	def __init__(self, wsheet_object):
		self._sheet = wsheet_object
		super(CheckmWorkflow, self).__init__(wsheet_object)
		options = [
			{"name": "input_file", "type": "infile", "format": "sequence.fasta"},
			{"name": "input_format", "type": "string", "default": "assembly"},
			{"name": "method", "type": "string", "default": "marker set"},
			{"name": "main_id", "type": "string"},
			{"name": "update_info", "type": "string"}
		]
		self.add_option(options)
		self.revise_infiles()
		self.set_options(self._sheet.options())

	def run(self):
		self.run_checkm()
		super(CheckmWorkflow, self).run()

	def run_checkm(self):
		self.checkm = self.add_tool("tool_lab.checkm")
		options = {
			"input_file": self.option('input_file'),
			"input_format": self.option('input_format'),
			"method": self.option('method'),
		}
		self.checkm.set_options(options)
		self.checkm.on('end', self.set_output)
		self.checkm.run()

	def set_output(self):	
		for file in os.listdir(self.checkm.output_dir):
			os.link(os.path.join(self.checkm.output_dir, file), os.path.join(self.output_dir, file))
		self.set_db()
		
	
	def set_db(self):
		for file in os.listdir(self.checkm.output_dir):
			self.logger.info("开始导表")
			api_checkm = self.api.api("tool_lab.checkm")
			api_checkm.add_Checkm(os.path.join(self.output_dir + "/" + file),main_id = self.option('main_id'))
			self.logger.info("导表结束")
		self.end()

	def end(self):
		result_dir = self.add_upload_dir(self.checkm.output_dir)
		super(CheckmWorkflow, self).end()