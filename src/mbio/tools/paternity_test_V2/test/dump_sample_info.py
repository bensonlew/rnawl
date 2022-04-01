# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'
import datetime
import pymongo
import codecs

analysistype = {"NIPT文库构建实验流程":"nipt", 
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
for k,v in analysistype.items():
	new = k.encode("utf-8").decode("utf-8")
	analysistype.pop(k)
	analysistype[new] = v

client = pymongo.MongoClient('mongodb://192.168.10.189:27017/')
db = client['tsanger_paternity_test_v2']
collection = db['sg_sample_info']

with codecs.open("/mnt/ilustre/users/sanger-dev/workspace/20171109/Single_file_check_test/FileCheckV4/20171019.interval.csv","r","utf-8") as INTER:
	title = INTER.readline()
	for line in INTER:
		linetemp = line.split(",")
		if linetemp[7] == "":
			pattern = linetemp[6]
		else:
			pattern = linetemp[6] + "-" + linetemp[7]
		print(pattern)
		insert_table = {
			'number':linetemp[0],
			'sample_accept_time':linetemp[1],
			'sample_name':linetemp[2],
			'library_type':linetemp[3],
			'analysis_type':analysistype[linetemp[3]],
			'product_type':linetemp[4],
			'sample_number':linetemp[5],
			'case_name':linetemp[6],
			'sample_id':linetemp[7],
			'sample_type':linetemp[8],
			'library_name':linetemp[9],
			's_sequence_num':linetemp[10],
			'index':linetemp[11],
			'sequence_size':linetemp[12],
			'board_batch':'test',
			'pattern':pattern,
		}
		collection.update({'pattern':pattern},{'$set':insert_table},upsert=True)