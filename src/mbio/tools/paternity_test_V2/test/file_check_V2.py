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
from biocluster.core.exceptions import OptionError
from optparse import OptionParser
from collections import Counter
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import MySQLdb as mdb
import numpy as np
import pdb
import codecs
import chardet
import sys
import os
import re
import shutil

class FileCheckAgent(Agent):
	def __init__(self,parent):
		super (FileCheckAgent,self).__init__(parent)
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


class FileCheckTool(Tool):
	def __init__(self, config):
		super(FileCheckTool,self).__init__(config)
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
		super(FileCheckTool, self).run()

		linestrip = self.linestrip
		DicConvertUTF8 = self.DicConvertUTF8
		Checkencoding = self.Checkencoding
		analysistype = self.analysistype
		sampletype = self.sampletype
		libraryType = self.libraryType

		DicConvertUTF8(sampletype)
		DicConvertUTF8(libraryType)
		DicConvertUTF8(analysistype)

		self.mergeupdown()
		self.Checkbarcode()
		self.NamePattern()
		self.Checkemergency()
		self.GenerateSplitTable()
		self.CheckSqlLIMS()
		self.end()

	def DicConvertUTF8(Dic):
		for k,v in Dic.items():
			new = k.encode("utf-8").decode("utf-8")
			Dic.pop(k)
			Dic[new] = v

	def Checkencoding(filepath):
		with open(filepath,"r") as f:
			title = f.readline()
			result = chardet.detect(title)
			if result['encoding'] != "utf-8":
				code = "gbk"
			else:
				code = "utf-8"

			return(code)

	def linestrip(line):
		line = line.strip()
		line = line.strip("\r")
		return(line)

	def GenerateSplitTable(self):
		out_title = ["sample_name","analysis_type","emergency","department","M reads","index","length"]
		IN = codecs.open(self.interval_add_emergency,"r","utf-8")
		OUT = codecs.open(self.split,"w","utf-8")
		IN_title = IN.readline()
		print(",".join(out_title),file=OUT)
		for line in IN:
			line = linestrip(line)
			linetemp = line.split(",")
			if linetemp[7] == "":
				sample_name = linetemp[6]
			else:
				sample_name = "-".join(linetemp[6:8])

			reformat = [sample_name,analysistype[linetemp[3]],linetemp[-1],"MED",linetemp[-4],linetemp[-3],linetemp[-2]]
			print(",".join(reformat),file=OUT)

		IN.close()
		OUT.close()

	def mergeupdown(self):
		###简单的用字典合并
		mergeDic = {}
		count = 1
		if(self.up != "empty"):
		#with codecs.open ("20170924.p.csv","r","gbk") as IN_UP:
			code_type = Checkencoding(self.up)
			with codecs.open (self.up,"r",code_type) as IN_UP: 
				title = linestrip(IN_UP.readline())
				titletemp = title.split(",")[0:14]
				titletemp[-1] = "备注"
				titletemp = ",".join(titletemp)
				mergeDic["title"] = titletemp
				for line in IN_UP:
					line = linestrip(line)
					linetemp = ",".join(line.split(",")[0:14])
					mergeDic[count] = linetemp
					count += 1

		#print ("up count is" + str(count))

		if(self.down != "empty"):
			#with codecs.open ("20170924.t.csv","r","gbk") as IN_DOWN:
			code_type = Checkencoding(self.down)
			#print ("down "+code_type)
			with codecs.open (self.down,"r",code_type) as IN_DOWN:
				IN_DOWN.next()
				IN_DOWN.next() ###前面三行没有作用
				IN_DOWN.next()
						
				for line in IN_DOWN:
					line = linestrip(line)
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

	def Checkbarcode(self):
		"""
			这一部分要处理简单的查重
			1. 首先检查barcode列是不是对的。
			2. 其次检查里面有没有重复
			3. 如果有重复，就要检查样品名是不是一致，如果一致就跳过。否则就要报错
		"""

		BarcodeDic = {}
		#with codecs.open("20171018.csv","r","utf-8") as IN:
		code_type = Checkencoding(self.interval)
		with codecs.open(self.interval,"r",code_type) as IN:
			title = IN.readline()
			for line in IN:
				line = linestrip(line)
				linetemp = line.split(",")
				if re.match(r'^[ATCG]+$',linetemp[-3]) == None:
					print("barcode consist of ATCG only,so barcode column may be error,check your file")
					sys.exit(1)
				if linetemp[-3] not in BarcodeDic:
					BarcodeDic[linetemp[-3]] = [linetemp[6]]
				else:
					BarcodeDic[linetemp[-3]].append(linetemp[6])
					samplename = list(set(BarcodeDic[linetemp[-3]]))
					if len(samplename) > 1:
						samplenameout = ",".join(samplename)
						print (samplenameout + " has same barcode " + linetemp[-3] + "but differ sample name")

		print("barcode check pass")

	def NamePattern(self):
		"""
			这一部分就是处理杂捕，多重之类的命名规则。
			字典放置在程序的最上层。考虑到以后很多东西都会重新命名
		"""
		error_count = 1
		code_type = Checkencoding(self.interval)
		#with codecs.open ("20171019.csv","r","utf-8") as IN:
		with codecs.open (self.interval,"r",code_type) as IN:
			title = IN.readline()
			for line in IN:
				line = linestrip(line)
				linetemp = line.split(",")

				#print(line)
				###检查样品类型是不是对的
				if sampletype.has_key(linetemp[8]):
					if re.search(r'^'+sampletype[linetemp[8]],linetemp[5]) == None:
						print (linetemp[5] + " do not have right sampletype " + linetemp[8] +" , check your files")
						error_count += 1
				else:
					continue

				###检查建库类型是不是对的
				if libraryType.has_key(linetemp[3]):
					if re.search(r'^'+libraryType[linetemp[3]],linetemp[6]) == None:
						print (linetemp[5] + " with differ librarytype " + linetemp[3])
						error_count += 1
				else:
					print (linetemp[3] + " do not have right librarytype, check your files")
					error_count += 1

			if error_count != 1:
				sys.exit(1)

	def Checkemergency(self):
		"""
			从样本批次表里面获取加急信息。
			1. 如果样本里没有批次信息，报错
			2. 批次信息数目只能 >= 样品数目。如果不是，就需要把这个表格生成
		"""
		Dic_PICI={}   ###保存批次表所有结果，在循环sample表的时候不断删除，最后得到批次表中有，sample表中没有的结果
		Dic_PICI_pass = {} ###保存sample表中有，batch表中也有的最终结果
		Dic_PICI_NOTIN_sample = {} ###保存sample表中有，而批次表中没有的结果，这个会报错

		interval_code = Checkencoding(self.interval)
		pici_code = Checkencoding(self.batch)
		INTERVAL = codecs.open(self.interval,"r",interval_code) 
		PICI = codecs.open(self.batch,"r",pici_code)

		PICI_title = linestrip(PICI.readline())
		INTERVAL_title = linestrip(INTERVAL.readline())
		for line in PICI:
			line = linestrip(line)
			linetemp = line.split(",")
			if Dic_PICI.has_key(linetemp[0]):
				print (linetemp[0] + " is replicated in batch table!!!check your files!!!")
			else:
				Dic_PICI[linetemp[0]] = line
				

		###开始循环interval表
		for line in INTERVAL:
			line = linestrip(line)
			linetemp = line.split(",")
			if linetemp[7] == "":
				pattern = linetemp[6]
			else:
				pattern = "-".join(linetemp[6:8])

			if Dic_PICI.has_key(pattern):
				linetemp[-1] = "" ##有时候拆分表会备注一些奇奇怪怪的东西。
				linetemp[-1] = (Dic_PICI[pattern].split(","))[-1]
				Dic_PICI_pass[pattern] = ",".join(linetemp)
				Dic_PICI.pop(pattern)
			else:
				Dic_PICI_NOTIN_sample[pattern] = line

		###开始判断结果
		if len(Dic_PICI_NOTIN_sample) > 0:
			print ("error:there are samples not in the batch table!!!")
			for k,v in Dic_PICI_NOTIN_sample.items():
				print (v)
			sys.exit(1)

		if len(Dic_PICI) > 0:
			print ("warning:there are samples in the batch table not sequenced!!!")

		INTERVAL.close()
		PICI.close()

		###输出结果
		with codecs.open(self.interval_add_emergency,"w","utf-8") as OUT:
			print (INTERVAL_title,file=OUT)
			for k in sorted(Dic_PICI_pass.keys()):
				print (Dic_PICI_pass[k],file=OUT)

	def CheckSqlLocal(self):

		#con = mdb.connect(dbconfigure["host"], dbconfigure["name"], dbconfigure["passwd"], dbconfigure["database"], init_command="set names utf8")
		con = mdb.connect("192.168.10.113","fengbo.zeng","die","gada", init_command="set names utf8") ###暂时先用本地数据库替代
		cur = con.cursor()
		cur.execute("select name from samples")
		rows = cur.fetchall()

		nameDic = {}
		count = 1
		for each in rows:
			nameDic[each[0]] = 0

		code_type = Checkencoding(self.interval_add_emergency)
		with codecs.open(self.interval_add_emergency,"r",code_type) as IN:
			title = IN.readline()
			for line in IN:
				line = linestrip(line)
				linetemp = line.split(",")
				if re.search(r'^WQ',linetemp[6]) != None:
					linename = "-".join(linetemp[6:8])
				elif re.search(r'^WS',linetemp[6]) != None:
					linename = linetemp[6]
				else:
					linename = None

				if linename != None and not nameDic.has_key(linename):
						print (linename + " is replicated with database check your files!")
						count += 1

		if count == 1:
			print("all name in " + self.interval + " check passed")
		else:
			sys.exit(1)

	def CheckSqlLIMS(self):
		"""
			这个模块主要执行3个功能
			1. 检查WQ的命名问题，必须包含F/M/S
			2. 检查亲子样品是否都有case
			3. 检查产筛样品是否都有case
		"""
		DicWQ = {}
		DicWS = {}
		WQ_sample_name = []
		WQ_error_count = 0
		with codecs.open("20171019.split.csv","r","utf-8") as IN:
		#with codecs.open(self.interval_add_emergency,"r","utf-8") as IN:
			for line in IN:
				line = linestrip(line)
				linetemp = line.split(",")
				if linetemp[0].startswith("WQ"):
					if re.search(r'-S|-M|-F',linetemp[0]) == None: ## -s -f -m
						WQ_error_count += 1
					case_name = linetemp[0].split("-")[0]
					DicWQ[case_name] = 1
					WQ_sample_name.append(linetemp[0])
				if linetemp[0].startswith("WS"):
					DicWS[linetemp[0]] = 1

		###链接数据库
		con = mdb.connect(dbconfigure["host"], dbconfigure["name"], dbconfigure["passwd"], dbconfigure["database"], charset = 'utf8')
		cur = con.cursor()

		if len(DicWQ) > 0:
			WQ_sample = "','".join(DicWQ.keys())
			sql = "SELECT i.note, o.create_date, o.name, i.patient_name, st.name, i.code, i.parented, o.sjrq, " \
					"i.accept_date, o.isurgent, o.gestational_weeks, t.name, t.principal FROM sample_info i " \
					"LEFT JOIN sample_order o ON i.sample_order = o.id LEFT JOIN primary_task t ON t.id = o.advance " \
					"LEFT JOIN dic_sample_type st ON st.id = i.sample_kind WHERE o.id IS NOT NULL AND i.note LIKE 'WQ%'" \
					"and i.note in ('{}') ORDER BY o.create_date DESC".format(WQ_sample)
			cur.execute(sql)
			case_infor = cur.fetchall()
			###循环核对样品名
			Dic_case = {}
			circle_error = 0
			for each in case_infor:
				eachtemp = "{}-{}".format(each[0],each[6])
				Dic_case[eachtemp] = 1
				###做一些特例,这部分的原则是，尽量让LIMS录入改好。以后碰到什么问题再修改对应的规则
				if len(each[6]) == 1:  ##如果长度只有1，只会出现两种情况，M或者F，把F1和F-1或者M1和M-1加入考虑
					eachtemp1 = "{}-{}1".format(each[0],each[6])
					eachtemp2 = "{}-{}-1".format(each[0],each[6])
					Dic_case[eachtemp1] = 1
					Dic_case[eachtemp2] = 1
				if len(each[6]) == 2:  ##如果长度是2的话，那么一般情况下是F1和M1，那么把F-1和M-1加入考虑
					eachtemp3 = "{}-{}-{}".format(each[0],each[6][0],each[6][1])
					Dic_case[eachtemp3] = 1
			for each in WQ_sample_name:
				if re.search(r'-S',each) == None:
					if re.search(r'-T',each) != None:
						each = each.split("-T")[0]
					if not Dic_case.has_key(each):
						print (each + " do not match database,check your files")
						circle_error += 1
			if circle_error >0:
				sys.exit(1)

