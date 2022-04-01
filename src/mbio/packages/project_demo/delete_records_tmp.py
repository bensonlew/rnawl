# -*- coding: utf-8 -*-
# __author__ = 'gdq'
import time
from biocluster.config import Config
import json
from bson import CodecOptions, SON

project_type='ref_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
collections = db.collection_names()
wf = open('sub2name.txt','w')
for coll in collections:
    if coll.endswith('_detail'):
        continue
    doc_dicts = db[coll].find().sort([('_id', -1)])
    doc_dict = list()
    try:
        doc_dict = list(doc_dicts)[0]
    except:
        pass
    print(doc_dict)
    if doc_dict:
        if 'params' in doc_dict:
            params_str = doc_dict['params']
            print(params_str)
            if params_str and params_str != "null" and len(params_str) > 1:
                if isinstance(params_str, dict):
                    if 'submit_location' in params_str:
                        wf.write('\"%s\":\"%s\",\n'%(doc_dict['params']['submit_location'],coll))
                else:
                    try:
                        params_dict = json.loads(params_str)
                        if 'submit_location' in params_dict:
                            wf.write('\"%s\":\"%s\",\n' % (params_dict['submit_location'], coll))
                    except:
                        pass
wf.close()

