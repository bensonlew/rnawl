# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'

import datetime
import pymongo
import codecs

def dump_error(self,status):
	with codecs.open("/mnt/ilustre/users/sanger-dev/workspace/20171107/Single_file_check_test/FileCheckV4/20171019.error.log","r","utf-8") as ERROR:
		errorlog = ERROR.readlines()
		errorlog = "".join(errorlog)

	client = pymongo.MongoClient('mongodb://192.168.10.189:27017/')
	db = client['tsanger_paternity_test_v2']

	insert_data = {
		'created_ts':str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
		'status':status,
		'desc':errorlog,
		'member_id':'4545',
	}
	collection = db['sg_file_check']
	run_id = collection.insert_one(insert_data).inserted_id
	print "add sg_file_check!"
	return run_id
