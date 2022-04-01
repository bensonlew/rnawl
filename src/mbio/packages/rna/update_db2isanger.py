# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import re
import sys
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
import json
import unittest
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from pymongo import MongoClient
import argparse


class UpdateDb2isanger(object):
    def __init__(self, project_type):
        super(UpdateDb2isanger, self).__init__()
        self.db = Config().get_mongo_client(mtype=project_type, db_version=1)[Config().get_mongo_dbname(project_type, db_version=1)]
        self._project_type = project_type
        self.sanger_ip = "10.11.1.102"
        #self._db_name = Config().MONGODB + '_ref_rna'

    def collection_backup(self, coll_name):
        '''
        备份collection到文件
        '''
        coll = self.db[coll_name]
        record_list = list()
        for record_dict in coll.find():
            record_dict.pop("_id")
            record_list.append(record_dict)
        with open(self._project_type + "." + coll_name, 'w') as f:
            f.write(json.dumps(record_list, indent=4, ensure_ascii=False))


    def updata_to_isanger(self, coll_name, json_file=None, t="replace"):
        "同步测试机collection到线上 新库"
        # isanger_client = MongoClient("mongodb://rna:y6a3t5n1w0y7@10.100.1.10/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
        # i_db = isanger_client.sanger_ref_rna_v2
        # i_collection = i_db[coll_name]
        if coll_name not in ["sg_table_relation", "sg_software_database", "software_database"]:
            raise Exception("{} can not delete".format(coll_name))
        collection = self.db[coll_name]
        
        try:
            if t == "replace":
                collection.delete_many({})
            with open(json_file, "r") as f:
                record_list = json.load(f)
            for record_dict in record_list:
                collection.insert_one(SON(record_dict))
        except Exception as e:
            print(e) 
            raise Exception("更新失败")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception("imput collection name")

    parser = argparse.ArgumentParser(description='for update by table\n ')
    parser.add_argument('-p', type=str, default="ref_rna_v2", help="project type.")
    parser.add_argument('-m', type=str, default=None, help='update/backup.')
    parser.add_argument('-j', type=str, default=None, help='json file.')
    parser.add_argument('-c', type=str, default=None, help='collection name.')
    

    args = parser.parse_args()
    os.environ["current_mode"]="workflow"
    os.environ["NTM_PORT"]="7322"
    os.environ["WFM_PORT"]="7321"
    api = UpdateDb2isanger(args.p)
    if args.m == "update":
        api.updata_to_isanger(args.c, json_file=args.j)
    elif args.m == "backup":
        api.collection_backup(coll_name=args.c)
