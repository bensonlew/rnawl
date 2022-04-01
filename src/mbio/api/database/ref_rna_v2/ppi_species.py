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
from api_base import ApiBase


class PpiSpecies(ApiBase):
    def __init__(self, bind_object):
        super(PpiSpecies, self).__init__(bind_object)
        project_type = "ref_rna_v2"
        self.ref_rna_db = Config().get_mongo_client(mtype=project_type, ref=True, db_version=1)[Config().get_mongo_dbname(project_type, ref=True, db_version=1)]


    def add_species_db(self, spe_path):
        '''
        添加数据库物名称信息
        '''
        collection = self.ref_rna_db["ppi_species"]
        data_list = []
        with open(spe_path, 'r') as f:
            for line in f.readlines():
                cols = line.strip().split("\t")
                data = {
                    "species_class": cols[0],
                    "species_id": cols[1],
                    "species_name": cols[2]
                }
                data_list.append(SON(data))
        if data_list:
            try:
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入ppi species信息失败")





if __name__ == '__main__':
    os.environ["current_mode"]="workflow"
    os.environ["NTM_PORT"]="7322"
    os.environ["WFM_PORT"]="7321"
    PpiSpecies(None).add_species_db(sys.argv[1])
