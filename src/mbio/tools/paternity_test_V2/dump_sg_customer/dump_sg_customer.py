#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
作者：kefei.huang
时间：20171205

"""
from __future__ import print_function
import codecs
import pymongo
import datetime
from biocluster.config import Config

def linestrip(line):
	line = line.strip()
	line = line.strip("\r")
	return(line)

def dump_sg_customer(self,MAINID):
	PICI_dic = {}
	SG = codecs.open("20171019.SG.csv","w","utf-8")
	with codecs.open("/mnt/ilustre/users/sanger-dev/sg-users/kefei.huang/file_check_test/20171019_pici.csv","r","gbk") as PICI:
		title = PICI.readline()
		for line in PICI:
			line = linestrip(line)
			linetemp = line.split(",")
			PICI_dic[linetemp[0]] = line

	db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
	collection = db['sg_customer']
	count = 1;
	with codecs.open("/mnt/ilustre/users/sanger-dev/workspace/20171205/Single_file_check_test/FileCheckV4/20171019.csv","r","utf-8") as INTER:
		title = INTER.readline()
		for line in INTER:
			line = linestrip(line)
			linetemp = line.split(",")
			if linetemp[7] == "":
				sample = linetemp[6]
			else:
				sample = "-".join(linetemp[6:8])
			##重构输出数据
			PICI_result = PICI_dic[sample].split(",")
			#tempout = [[sample,PICI_result[1:3],linetemp[3],linetemp[9:14]]]
			tempout = [sample]
			tempout.extend(PICI_result[1:3])
			tempout.append(linetemp[3])
			tempout.extend(linetemp[9:14])
			tempout.append(PICI_result[3])
			tempout = ",".join(list(tempout))
			print(tempout,file=SG)

			insert_table = {
				"number" : count,
				"sample_id" : tempout[0],
				"extract_batch" : tempout[1],
				"library_batch" : tempout[2],
				"library_type" : tempout[3],
				"library_name" : tempout[4],
				"index" : tempout[6],
				"s_sequence_num" : tempout[5],
				"sequence_size" : tempout[7],
				"urgence" : tempout[9],
				"use_type" : tempout[8],
				"upload_time" : datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
				"board_batch" : "test",
				"check_id" : ObjectId(MAINID)
			}
			count += 1
			collection.update({'board_batch':self.flowcell},{'$set':insert_table},upsert = True)
	SG.close()

