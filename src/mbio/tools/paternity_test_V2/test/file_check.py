## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "kefei.huang"
#last_modify:20171011

import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import pandas as pd
import os
import re


"""
程序对输入的两个文件：线上拆分表，线下拆分表进行检查。程序逻辑如下
1. 合并两个文件
2. 与SQL数据库检查，有没有重复名称，如果有跳出并报错
3. 检查index是否有重复，如果有重复且不是同一个编号，跳出报错
4. 检查建库类型是否一致
"""

class FileCheckAgent(Agent):
	def __init__(self,parent):
		super (FileCheckAgent,self).__init__(parent)
		self.logger.info("run")
		options = [
			{"name":"up","type":"string"},
			{"name":"down","type":"string"}
			#{"name":"up","type":"string"},
			#{"name":"down","type":"string"}
		]
		self.add_option(options)
		self.logger.info("options finished")

	def check_options(self):
		if not self.option("up"):
			raise OptionError("请输入线上表格")
		if not self.option("down"):
			raise OptionError("请输入线下表格")
		return True

	def set_resource(self):
		self._cpu = 1
		self._memory = "1G"

class FileCheckTool(Tool):
	def __init__(self, config):
		super(FileCheckTool,self).__init__(config)

	def run(self):
		super(FileCheckTool, self).run()
		self.up_down_check()
		self.end()

	#def run_check(self):
		#print("check finished")


	def up_down_check(self):
		#up_pd = pd.read_csv("20170924.p.csv",encoding="gbk",index_col=False)
		#down_pd = pd.read_csv("20170924.t.csv",encoding="gbk",header=2,index_col=False)

up_pd = pd.read_csv(self.option("up"),encoding="GBK")
down_pd = pd.read_csv(self.option("down"),encoding="gbk",header=2)
down_pd = down_pd.iloc[:,0:11]
down_pd = down_pd.dropna(axis=0,how='all') ###保留11列，去除其他NA行
down_temp = down_pd.iloc[:,[5,5,6]]
down_temp.columns = ['序号','收样日期','样本姓名']
down_pd = pd.concat([down_temp,down_pd],axis=1)
up_pd = up_pd.iloc[:,0:14]
up_pd = up_pd.dropna(axis=0,how='all')
pd.concat([up_pd,down_pd],axis=0,join='outer')
