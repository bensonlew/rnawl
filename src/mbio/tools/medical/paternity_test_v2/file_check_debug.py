#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
作者：kefei.huang
时间：20171023

这个程序是为了合并三个数据表格而写的

1. 先把两个表格合并。线上表格是必须会有的，线下可能没有
2。随后对数据进行检查
3. 再调用批次表，得到加急非加急
4. 再调用数据库进行查重
5. 最后导表

"""
from __future__ import print_function
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import re
import codecs
import chardet


class FileCheckDebugAgent(Agent):
	def __init__(self,parent):
		super (FileCheckDebugAgent,self).__init__(parent)
		self.logger.info("run")
		options = [
			{"name":"up","type":"string"},
			{"name":"down","type":"string"},
			{"name":"date","type":"string"},
			{"name":"batch","type":"string"}
		]
		self.add_option(options)
		self.logger.info("options finished")

	def check_options(self):
		if not self.option("up"):
			raise OptionError("请输入线上表格")
		if not self.option("down"):
			raise OptionError("请输入线下表格")
		if not self.option("date"):
			raise OptionError("请输入结果文件名称，以时间为单位，格式请参照20171023")
		if not self.option("batch"):
			raise OptionError("请输入实验批次表")
		return True

	def set_resource(self):
		self._cpu = 1
		self._memory = "1G"

	def end(self):
		super(FileCheckDebugAgent,self).end()


class FileCheckDebugTool(Tool):
	def __init__(self, config):
		super(FileCheckDebugTool,self).__init__(config)
		self.up = self.option("up")
		self.down = self.option("down")
		self.date = self.option("date")
		self.batch = self.option("batch")
		self.interval = self.date + ".csv"
		self.interval_add_emergency = self.date + ".interval.csv"
		self.split = self.date + ".split.csv"

		self.dbconfigure={'host':'172.16.101.202','name':'mjlimsbak','passwd':'Q6k9pU9ZN5FTow','database':'mjlims'}
		self.libraryType = {"NIPT文库构建实验流程":"WS", 
									"WQ gDNA建库实验流程_多重":"WQ", 
									"WQ cfDNA建库实验流程_杂捕":"WQ",
									"FA gDNA建库实验流程_杂捕":"FA",
									"FA ctDNA建库实验流程_杂捕":"FA",
									"FA gDNA建库实验流程":"FA",
									"FA cfDNA建库实验流程_杂捕":"FA",
									"YCZL 建库实验流程":"YCZL",
									"QA 建库实验流程":"QA",
									"CT 建库实验流程":"CT",
									"HLY 建库实验流程":"HLY",
									"RXA 建库实验流程":"RXA",
									"snp 建库实验流程":"SNP",
									"XA RNA建库实验流程":"XA",
									"XA DNA建库实验流程":"XA",
									"甲基化 建库实验流程":"METH",
									"small RNA建库实验流程":"small"}
		self.sampletype = {"全血": "QX", 
									"血浆": "XJ", 
									"蜡块": "SL", 
									"石蜡": "SL", 
									"石蜡切片": "SL", 
									"石蜡切片(白片)": "SL", 
									"胸腹水" : "XS", 
									"手术标本": "XZ", 
									"穿刺标本": "XZ", 
									"穿刺样本": "XZ",
									"组织标本": "XZ", 
									"新鲜组织": "XZ",
									"蜡卷": "SL",
									"精斑": "JB",
									"亲子父本全血":"QQ",
									"指甲": "ZJ"}
		self.analysistype = {"NIPT文库构建实验流程":"nipt", 
									"WQ gDNA建库实验流程_多重":"dcpt", 
									"WQ cfDNA建库实验流程_杂捕":"pt",
									"FA gDNA建库实验流程_杂捕":"ctdna",
									"FA ctDNA建库实验流程_杂捕":"ctdna",
									"FA gDNA建库实验流程":"ctdna",
									"FA cfDNA建库实验流程_杂捕":"ctdna",
									"YCZL 建库实验流程":"genetic_tumor",
									"QA 建库实验流程":"QA_str",
									"CT 建库实验流程":"ct_str",
									"HLY 建库实验流程":"dc_str",
									"XA RNA建库实验流程":"blood_cancer",
									"XA DNA建库实验流程":"blood_cancer",
									"甲基化 建库实验流程":"",
									"RXA 建库实验流程":"",
									"snp 建库实验流程":"",
									"small RNA建库实验流程":""}

	def run(self):
		super(FileCheckDebugTool, self).run()
		self.logger.info("start run")
		self.DicConvertUTF8()
		self.mergeupdown()
		#self.Checkbarcode()
		#self.NamePattern()
		#self.Checkemergency()
		#self.GenerateSplitTable()
		#self.CheckSqlLIMS()
		self.end()

	def DicConvertUTF8(self):
		for k,v in self.sampletype.items():
			new = k.encode("utf-8").decode("utf-8")
			self.sampletype.pop(k)
			self.sampletype[new] = v

		for k,v in self.libraryType.items():
			new = k.encode("utf-8").decode("utf-8")
			self.libraryType.pop(k)
			self.libraryType[new] = v

		for k,v in self.analysistype.items():
			new = k.encode("utf-8").decode("utf-8")
			self.analysistype.pop(k)
			self.analysistype[new] = v

	def Checkencoding(self,filepath):
		with open(filepath,"r") as f:
			title = f.readline()
			result = chardet.detect(title)
			if result['encoding'] != "utf-8":
				code = "gbk"
			else:
				code = "utf-8"

			return(code)

	def linestrip(self,line):
		line = line.strip()
		line = line.strip("\r")
		return(line)

	def mergeupdown(self):
		###简单的用字典合并
		mergeDic = {}
		count = 1
		if(self.up != "empty"):
		#with codecs.open ("20170924.p.csv","r","gbk") as IN_UP:
			code_type = self.Checkencoding(self.up)
			with codecs.open (self.up,"r",code_type) as IN_UP: 
				title = self.linestrip(IN_UP.readline())
				titletemp = title.split(",")[0:14]
				titletemp[-1] = "备注"
				titletemp = ",".join(titletemp)
				mergeDic["title"] = titletemp
				for line in IN_UP:
					line = self.linestrip(line)
					linetemp = ",".join(line.split(",")[0:14])
					mergeDic[count] = linetemp
					count += 1

		#print ("up count is" + str(count))

		if(self.down != "empty"):
			#with codecs.open ("20170924.t.csv","r","gbk") as IN_DOWN:
			code_type = self.Checkencoding(self.down)
			#print ("down "+code_type)
			with codecs.open (self.down,"r",code_type) as IN_DOWN:
				IN_DOWN.next()
				IN_DOWN.next() ###前面三行没有作用
				IN_DOWN.next()
						
				for line in IN_DOWN:
					line = self.linestrip(line)
					linetemp = line.split(",")[0:11]
					if linetemp[0] == "" and linetemp[1] == "":
						break

					pushtemp = ["","",""]
					pushtemp.extend(linetemp)
					mergeDic[count] = ",".join(pushtemp)
					count += 1
					#print (count)

		#print ("down count is" + str(count))

		with codecs.open(self.interval,"w","utf-8") as out:
			print(mergeDic["title"],file=out)
			for k in sorted(mergeDic.keys()):
				if k == "title" or mergeDic[k] == None:
					continue
				else:
					print(mergeDic[k],file = out)
