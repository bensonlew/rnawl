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
from mbio.api.database.denovo_rna_v2.api_base import ApiBase


class TaxonDb(ApiBase):
    def __init__(self, bind_object):
        super(TaxonDb, self).__init__(bind_object)
        #self._db_name = Config().MONGODB + '_ref_rna'

    def add_tax_db(self, name_dmp):
        '''
        添加数据库物名称信息
        '''
        collection = self.db["sg_taxon_db"]
        data_list = []
        with open(name_dmp, 'r') as f:
            for line in f.readlines():
                cols = line.split("\t")
                if cols[6] == "misspelling":
                    pass
                else:
                    species = {
                        "taxon_id": cols[0],
                        "species_name": cols[2],
                        "name_type": cols[6]
                    }
                    data_list.append(SON(species))
        if data_list:
            try:
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入taxon name信息失败" , code="52002701")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        self.bind_object.set_error("python genome_db.py name.dmp", code="52002702")
    TaxonDb(None).add_tax_db(sys.argv[1])
