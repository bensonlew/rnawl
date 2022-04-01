#!#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kefei.huang
# last modified 20171222 by kefei.huang
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import os

class SampleSimilarityModule(Module):
	def __init__(self, work_id):
		super(SampleSimilarityModule, self).__init__(work_id)
		options = [
			{"name": "sample_id", "type": "string"}  # 输入F/M/S的样本ID
		]
		self.add_option(options)
		self.run_tools = []

	def check_options(self):
		if not self.option("sample_id"):
			raise OptionError("请输入样品名称")
		##检查文件名是不是都存在
		self.library_path = Config().SOFTWARE_DIR + "/" + "database/human/pt_ref/tab_all"		
		all_samples = self.option("sample_id").split(",")
		for m in all_samples:
			abs_sample = "{}/{}.tab".format(self.library_path,m)
			if not os.path.exists(abs_sample):
				raise OptionError("{}文件不存在,请确认".format(abs_sample))
		return True

	def finish_update(self, event):
		step = getattr(self.step, event['data'])
		step.finish()
		self.step.update()

	def similarity_run(self):
		n = 0
		all_samples = self.option("sample_id").split(",")
		for m in all_samples:
			single_run = self.add_tool("medical.paternity_test_v2.sample_similarity")
			self.step.add_steps('sample_similarity_{}'.format(n))
			single_run.set_options({
				"sample_id":m				
			})
			step = getattr(self.step, 'sample_similarity_{}'.format(n))
			step.start()
			single_run.on('end', self.finish_update, 'sample_similarity_{}'.format(n))
			self.run_tools.append(single_run)
			n += 1
			self.on_rely(self.run_tools, self.end)
		for tools in self.run_tools:
		    tools.run()

	def run(self):
		super(SampleSimilarityModule, self).run()
		self.similarity_run()

	def end(self):
		super(SampleSimilarityModule, self).end()

