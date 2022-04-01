# -*- coding: utf-8 -*-
# __author__ :zhouxuan

import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re

class DataSplitAgent(Agent):

	"""
		项目：亲子鉴定
		功能：对医学检验所的测序数据进行拆分，区分各个数据(WQ对应亲子鉴定，WS对应产前筛查等)
		author：zhouxuan 2017.02.21
		last modify: 2017.04.21
		version: v1.0

	"""
	def __init__(self, parent):
		super(DataSplitAgent, self).__init__(parent)
		options = [
			{"name": "message_table", "type": "infile", "format": "paternity_test.tab"},
			{"name": "data_dir", "type": "string"},
			{"name": "ws_single", "type": "string"},
		]
		self.add_option(options)
		self.step.add_steps("data_split")
		self.on('start', self.stepstart)
		self.on('end', self.stepfinish)

	def stepstart(self):
		self.step.data_split.start()
		self.step.update()

	def stepfinish(self):
		self.step.data_split.finish()
		self.step.update()

	def check_options(self):
		"""
		重写参数检测函数
		:return:
		"""
		if not self.option('message_table'):
			raise OptionError("必须输入样本信息表")
		if not self.option('data_dir'):
			raise OptionError("必须提供样本序列文件夹")
		if not self.option('ws_single'):
			raise OptionError("必须提供是否只有ws单端序列")
		return True

	def set_resource(self):  # 后续需要测试确认
		"""
		设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
		:return:
		"""
		self._cpu = 50
		self._memory = '100G'

	def end(self):
		result_dir = self.add_upload_dir(self.output_dir)
		result_dir.add_relpath_rules([
			[".", "", "结果输出目录"]
		])
		result_dir.add_regexp_rules([
			["", "", ""],
		])
		super(DataSplitAgent, self).end()

class DataSplitTool(Tool):
	def __init__(self, config):
		super(DataSplitTool, self).__init__(config)
		self._version = "v1.0"
		self.script_path = "bioinfo/medical/bcl2fastq-2.17/bin/bcl2fastq"
		self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
		self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')

	def run(self):
		"""
		运行
		:return:
		"""
		super(DataSplitTool, self).run()
		self.run_ds()
		self.set_output()
		self.end()

	def run_ds(self):
		"""
		处理输入文件并完成数据拆分
		:return:
		"""
		message_path = self.option('message_table').prop['path']
		new_message_table = os.path.join(self.work_dir, "new_message_table")
		with open(new_message_table, "a") as w:
			w.write("[Data],,,,,,," + "\n" +
			        "Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description" + "\n")
		with open(message_path, "r") as r:
			content = r.readlines()
			for line in content:
				line = line.split("\t")
				with open(new_message_table, "a") as w:
					lines = "Sample_" + line[3] + "," + line[3] + ",,,," + line[8] + "," + line[4] + "," + "\n"
					w.write(lines)

		index = self.option('data_dir').split(":")[0]
		name = Config().get_netdata_config(index)
		old_data_dir = name[index + "_path"] + "/" + self.option('data_dir').split(":")[1]
		if not os.path.exists(old_data_dir):
			self.set_error("下机数据文件夹路径不正确，请设置正确的路径。")
			raise Exception("下机数据文件夹路径不正确，请设置正确的路径。")
		if self.option('ws_single') == 'false':
			cmd = "{} -i {}Data/Intensities/BaseCalls/ -o {} --sample-sheet {} --use-bases-mask  y76,i6n,y76 " \
			      "--ignore-missing-bcl -R {} -r 4 -w 4 -d 2 -p 10 --barcode-mismatches 0".\
				format(self.script_path,old_data_dir,self.output_dir,
			           new_message_table, old_data_dir)
		else:
			cmd = "{} -i {}Data/Intensities/BaseCalls/ -o {} --sample-sheet {} --use-bases-mask  y76,i6n " \
			      "--ignore-missing-bcl -R {} -r 4 -w 4 -d 2 -p 10 --barcode-mismatches 0". \
				format(self.script_path, old_data_dir, self.output_dir,
			           new_message_table, old_data_dir)
		self.logger.info("start data_split")
		command = self.add_command("ds_cmd", cmd).run()
		self.wait(command)
		if command.return_code == 0:
			self.logger.info("data_split done")
		else:
			self.set_error("data_split error")
			raise Exception("data_split error")

	def set_output(self):
		"""
		将数据拆分的结果，直接存放在同一个文件夹下，按照样本分成不同的小文件夹
		:return:
		"""
		m = re.match('(.*)_A(.*)/', self.option('data_dir'))
		# n = re.match('(.*)_A(.*)_.*/', self.option('data_dir'))
		if m:
			chip_name = m.group(2)[0:9]
		else:
			raise Exception('无法获取到相应的数据的芯片名'.format(self.option('data_dir')))
		self.logger.info("芯片名字:{}".format(chip_name))
		origin_html = self.output_dir + '/Reports/html/' + chip_name +'/all/all/all/laneBarcode.html'
		if os.path.exists(origin_html):
			sample = self.check_yield(origin_html, self.work_dir + '/yield.txt')
			if len(sample) > 0:
				self.set_error("以下样本yield为0：{}".format(sample))
				raise Exception('以下样本yield为0：{}'.format(sample))
			else:
				self.logger.info('拆分成功')
		else:
			raise Exception("拆分失败")

	@staticmethod
	def check_yield(html_file, txt_table):  # 检查是否存在空的样本
		'''
		拆分结束后检查样本是否为空
		:param html_file: 拆分生成的文件
		:param txt_table: html导出的TXT表格
		:return: 数据量为空的样本
		'''
		start = 'false'
		with open(html_file, 'r') as h, open(txt_table, 'w') as t:
			t.write('Lane\tProject\tSample\tBarcode sequence\tPF Clusters\t% of the lane\t% Perfect barcode\t'
			        '% One mismatch barcode\tYield (Mbases)\t% PF Clusters\t% >= Q30 bases\tMean Quality Score')
			for line in h:
				content = line.rstrip()
				if content == '<th>Mean Quality<br>Score</th>':
					start = 'true'
				else:
					if start == 'false':
						continue
					else:
						if content == '</tr>':
							t.write('\n')
						elif content == '<h2>Top Unknown Barcodes</h2>':
							break
						else:
							m = re.match('<td>(.*)</td>', content)
							if m:
								t.write(m.group(1) + '\t')
							else:
								pass
		problem_sample = []
		with open(txt_table, 'r') as r:
			for line in r:
				content = line.rstrip().split('\t')
				if content[8] == '0':
					problem_sample.append(content[2])
		return set(problem_sample)

