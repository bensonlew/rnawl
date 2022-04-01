# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'

import datetime
import pymongo
import codecs
client = pymongo.MongoClient('mongodb://192.168.10.189:27017/')
db = client['tsanger_paternity_test_v2']
collection = db['sg_batch_test']
with codecs.open("/mnt/ilustre/users/sanger-dev/sg-users/kefei.huang/file_check_test/20171019_pici.csv","r","gbk") as PICI:
	title = PICI.readline()
	for line in PICI:
		linetemp = line.split(",")
		if linetemp[0].startswith('WQ') or linetemp[0].startswith('WS'):
			insert_data = {
				'sample_id':linetemp[0],
				'extract_batch':linetemp[1],
				'library_batch':linetemp[2],
				'urgency_degree':linetemp[3],
				'board_batch':'test',
				'insert_time':datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
			}
			status = collection.update({"sample_id":linetemp[0]},{'$set':insert_data},upsert=True)
