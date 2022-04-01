# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'
import os
import re
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import pymongo
import time
import codecs

#client = pymongo.MongoClient('mongodb://192.168.10.189:27017/')
#db = client['tsanger_paternity_test']
#test = db['sg_pt_qc'].find_one({'sample_id':'WQ17072887-M-1'})
#if test != None:
def CheckMongoDupName(self):
	error_count = 0;
	client = pymongo.MongoClient('mongodb://192.168.10.187:27017/')
	db = client['sanger_paternity_test_ref']

	with codecs.open("/mnt/ilustre/users/sanger-dev/workspace/20171107/Single_file_check_test/FileCheckV4/20171019.error.log","r","utf-8") as ERROR:
		errorlog = ERROR.readlines()

	with codecs.open("20171019.split.csv","r","utf-8") as IN:
		title = IN.readline()
		for line in IN:
			linetemp = line.split(',')
			name = linetemp[0]
			name_return = db['sg_pt_qc'].find_one({'sample_id':name})
			if name_return != None:
				print (name + " is already used,check your files and change name")
				error_count += 1
	return(error_count)

	