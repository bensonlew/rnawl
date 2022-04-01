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


class TaxonDb(ApiBase):
    def __init__(self, bind_object):
        super(TaxonDb, self).__init__(bind_object)


    def add_tax_db(self, name_dmp):
        '''
        添加数据库物名称信息
        '''
        collection = self.db["ncbi_taxon_name"]
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

    def search_taxon(self, species):
        species = species.replace("_", " ")
        spe_coll = 'ncbi_taxon_name'
        collect = self.db[spe_coll]
        score = {"score":{"$meta":"textScore"}}
        sort_list = [("score", {"$meta":"textScore"})]
        condition = {"$text": {"$search": species}}
        # target = {"taxon_id": 1, "species_name": 1}
        results = collect.find(condition, score).sort(sort_list).limit(20)

        #.sort(score).limit(20)
        spe_list = [(d["species_name"], d["taxon_id"]) for d in results]
        spe_dict = dict(spe_list)
        if species in spe_dict:
            return {species: spe_dict[species]}
        else:
            return dict(spe_list[0:1])

    def search_taxon_by_id(self, taxon_id):
        spe_coll = 'ncbi_taxon_name'
        collect = self.db[spe_coll]
        condition = {"taxon_id": taxon_id, "name_type": "scientific name"}
        # target = {"taxon_id": 1, "species_name": 1}
        result = collect.find_one(condition)
        if result:
            return result["species_name"]
        else:
            return None


if __name__ == '__main__':
    # TaxonDb(None).add_tax_db(sys.argv[1])
    TaxonDb(None).search_taxon(sys.argv[1])
