# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from __future__ import print_function
import pymongo
import json
from biocluster.config import Config

table_json = dict()
params_json = dict()
sub2name_json = dict()
params_json['common'] = {'task_id': '', 'project_sn': '', 'submit_location': '', 'task_type': ''}
# denovo_rna_v2
db = Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')]
params_json['small_rna'] = dict()
table_json['small_rna'] = dict()
sub2name_json['small_rna'] = dict()
table_relation = db['sg_table_relation'].find_one({})['target']
main_tables = list()
for target in table_relation:
    if target[0] not in main_tables:
        main_tables.append(target[0])
    if target[2] and target[2] not in table_json['small_rna']:
        table_json['small_rna'][target[2]] = target[0]
coll_names = db.collection_names()
for collection in coll_names:
    if collection not in main_tables:
        continue
    try:
        result = db[collection].find_one({"status": "end"}, sort=[('_id', pymongo.DESCENDING)])
    except:
        continue
    if result and 'params' in result and result['params']:
        try:
            params = json.loads(result['params'])
        except:
            continue
        if isinstance(params, dict) and 'submit_location' in params:
            submit_location = params['submit_location'].encode("utf-8")
            results = db['sg_status'].find({"status": "end", "submit_location": submit_location})
            total_params = list()
            for result in results:
                params_temp = json.loads(result['params'])
                for para in params_temp:
                    if para not in total_params:
                        total_params.append(para)
            if submit_location not in sub2name_json['small_rna']:
                sub2name_json['small_rna'][submit_location] = ""
            params_json['small_rna'][submit_location] = dict()
            for key in total_params:
                if key in ['task_id', 'project_sn', 'submit_location', 'task_type', 'client']:
                    continue
                if isinstance(key, unicode):
                    key = key.encode("utf-8")
                if key not in params_json['small_rna'][submit_location]:
                    params_json['small_rna'][submit_location][key] = ''
with open('sublocation2params.json', 'w') as f:
    f.write(json.dumps(params_json, indent=4))
with open('mainid2table.json', 'w') as f:
    f.write(json.dumps(table_json, indent=4))
with open('sub2name.json', 'w') as f:
    f.write(json.dumps(sub2name_json, indent=4))