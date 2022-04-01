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


class UpdateDb2isanger(ApiBase):
    def __init__(self, bind_object):
        super(UpdateDb2isanger, self).__init__(bind_object)
        #self._db_name = Config().MONGODB + '_ref_rna'

    def collection_backup(self, coll_name)
        record_list = list()
        

    def updata_to_isanger(self, coll_name, t="replace"):
        "根据genome_id获取genome属性字典"
        isanger_client = MongoClient("mongodb://rna:y6a3t5n1w0y7@10.100.1.10/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
        i_db = isanger_client.sanger_ref_rna_v2
        i_collection = i_db[coll_name]
        if coll_name not in ["sg_table_relation", "sg_software_database"]:
            raise Exception("{} can not delete".format(coll_name))
        collection = self.db[coll_name]
        try:
            if t == "replace":
                with open(coll_name + "old.txt", 'w') as f:
                    for old_dict in i_collection.find():
                        f.write("{}\n".format(old_dict))
                i_collection.delete_many({})
            record_dicts = collection.find()
            for record_dict in record_dicts:
                i_collection.insert_one(SON(record_dict))
        except Exception as e:
            print e
            raise Exception("更新失败")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception("imput collection name")

    UpdateDb2isanger(None).updata_to_isanger(sys.argv[1])
