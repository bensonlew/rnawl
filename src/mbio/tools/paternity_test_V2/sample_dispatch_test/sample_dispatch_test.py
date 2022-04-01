#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
作者：kefei.huang
时间：20171115

使用datasplit的结果，根据特定的规则，对样品的投递顺序进行排序。
目前主要有两个workflow，亲子和产筛。
目前来看，将这两个workflow设定为并行状态，不排除以后根据特别的情况，做成workflow串行的状态。
所以目前，根据投递时间还有投递的资源，粗略的控制投递顺序。

"""
from __future__ import print_function
from biocluster.config import Config
import codecs
import pymongo
import re
###首先打开datasplit的结果,拿到需要运行的家系
###整体的运行顺序，以后要修改顺序和增删就在这里
ALL_report = {
	'pt-emergercy':[],
	'pt-hurry':[],
	'pt-normal':[],
	'report-emergency':[],
	'report-hurry':[],
	'report-normal':[],
	'nipt-emergency':[],
	'nipt-hurry':[],
	'nipt-normal':[],
	'other':[]
}

emergencytype = {
	"pt-特加急":2,
	"pt-加急":1,
	"dcpt-特加急":2,
	"dcpt-加急":1,
	"nipt-加急":1,
	"nipt-特加急":2
}

for k,v in emergencytype.items():
	new = k.encode("utf-8").decode("utf-8")
	emergencytype.pop(k)
	emergencytype[new] = v

caseDic = {}
StatuDic = {}
with codecs.open("20171019.split.csv","r","utf-8") as SPLIT:
	title = SPLIT.readline()
	for line in SPLIT:
		linetemp = line.split(",")
		if (linetemp[2] == ""):
			pattern = linetemp[1]
		else:
			pattern = "-".join(linetemp[1:3])
		if linetemp[0].startswith("WQ"):
			casetemp = linetemp[0].split("-")
			StatuDic[linetemp[0]] = pattern
			if caseDic.has_key(casetemp[0]):
				caseDic[casetemp[0]].append(linetemp[0])
			else:
				caseDic[casetemp[0]] = [linetemp[0]]
		elif linetemp[0].startswith("WS"):
			if emergencytype.has_key(pattern):
				if emergencytype[pattern] == 2:
					ALL_report['nipt-emergency'].append(linetemp[0])
				elif emergencytype[pattern] == 1:
					ALL_report['nipt-hurry'].append(linetemp[0])
			else:
				ALL_report['nipt-normal'].append(linetemp[0])
		else:
			ALL_report['other'].append(linetemp[0])

###分类
###key:dcpt:F,M，是第一步要运行的数据
###key:pt-特加急与加急: 需要加急运行的家系信息
###key:nipt-特加急与加急：需要加急运行的产筛样品
###key:pt: 非加急运行的家系信息
###key:nipt 非加急运行的产筛样品
###key:其他  没有写明的都是不再线上运行的样品


###开始计算需要运行的家系
###1. 首先把所有的PT样品的CASE收集起来。得到需要穷举的家系
###2. 刨除已经运行完的家系
###3. 把家系分为特急 + 加急 + 普通
###4. 把特急 + 加急 + 普通 分成 直出报告 与 先运行S + 报告

###暂时使用 sg_pt_ref_main表，使用正则表达式来搜寻
#client = pymongo.MongoClient('mongodb://192.168.10.187:27017/')
#db = client['tsanger_paternity_test_ref']
#collection = db['sg_pt_ref_main']
db = Config().get_mongo_client(mtype="pt", ref=True)[Config().get_mongo_dbname("pt", ref=True)]
collection = db['sg_pt_ref_main']

for k in caseDic.keys():
	#result = collection.find({'sample_id':'WQ17093590-S-1')
	result=collection.find({'sample_id':{'$regex':k}})
	if result.count() > 0:
		result = list(result)
		for temp in result:
			caseDic[k].append(temp['sample_id'])

###此时拿到的数据是所有的样品,开始循环生成所有的报告
report = []
for k in caseDic.keys():
	father = []
	mother = []
	son = []
	for each in caseDic[k]:
		if re.search("-S",each):
			son.append(each)
		elif re.search("-M",each):
			mother.append(each)
		elif re.search("-F",each):
			father.append(each)
		else:
			print (each + " is not one of father,mother and son, check files")
	if (len(father) >0 and len(mother) >0 and len(son)>0):
		for son_each in son:
			for mother_each in mother:
				for father_each in father:
					name = "_".join([son_each,mother_each,father_each])
					report.append(name)

###我们已经拿到了所有的报告。接下来对报告分类
report = list(set(report))

for each_report in report:
	report_type = []
	each_temp = each_report.split("_")
	fast_count = 0;
	no_count = 0;
	for temp in each_temp:
		if StatuDic.has_key(temp):
			report_type.append(StatuDic[temp])
			if emergencytype.has_key(StatuDic[temp]):
				if fast_count < emergencytype[StatuDic[temp]]:
					fast_count = emergencytype[StatuDic[temp]]
		else:
			report_type.append("no") ###不在这份列表里面
			no_count += 1;
	###再循环一次开始做处理
	if no_count == 3:
		continue
	if fast_count == 2:
		if report_type[0] != "no":
			ALL_report['pt-emergercy'].append(each_report)
		else:
			ALL_report['report-emergency'].append(each_report)
	elif fast_count == 1:
		if report_type[0] != "no":
			ALL_report['pt-hurry'].append(each_report)
		else:
			ALL_report['report-hurry'].append(each_report)
	else:
		if report_type[0] != "no":
			ALL_report['pt-normal'].append(each_report)
		else:
			ALL_report['report-normal'].append(each_report)


###输出结果
with codecs.open("20171019_report.csv","w","utf-8") as OUT:
	ALL_report_keys = sorted(ALL_report.keys(),reverse=True)
	for k in ALL_report_keys:
		print ("\n"+k,file=OUT)
		print ("\n".join(ALL_report[k]),file=OUT)




