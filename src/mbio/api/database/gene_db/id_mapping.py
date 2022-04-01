# !/usr/bin/python
# -*- coding: utf-8 -*-
import types
import re
import os
import json
import math
from collections import OrderedDict
import unittest
import datetime
import glob
from bson.objectid import ObjectId
from api_base import ApiBase
import csv
import sys
import inspect
# import pandas as pd
from pymongo import IndexModel, ASCENDING, DESCENDING, TEXT


class IdMapping(ApiBase):
    def __init__(self, bind_object):
        super(IdMapping, self).__init__(bind_object)

    def add_pipe(self, species, gb_acc, species_name=None, species_name2=None, gene_id_file=None, tran_id_file=None, ensemble_schema=None, uniprot_schema=None):
        """
        自动添加数据
        """
        logic_col = self.db["sgdb_logic"]
        genome_col = self.db["sgdb_genome"]
        logic_table = self.get_tables_by_main_record("sgdb_logic")
        field_list = [dic["field"] for dic in logic_table]

        data_list = []
        with open(ensemble_schema, "r") as f:
            for dic in csv.DictReader(f, delimiter='\t'):
                if dic['internalName_c'] == "xrefs":
                    if dic['internalName'] in field_list:
                        pass
                    else:
                        if dic["linkoutURL"].startswith("exturl"):
                            if dic["linkoutURL"].startswith("exturl|/"):
                                url = "exturl|/" + "/".join(dic["linkoutURL"].split("/")[2:])
                            else:
                                url = dic["linkoutURL"][7:]
                        else:
                            url = dic["linkoutURL"]

                        data = {
                            "uniprot_match" : "",
                            "displayName" : dic["displayName"],
                            "NCBI" : "",
                            "url" : url,
                            "selected" : "ENSEMBL",
                            "uniprot_use" : "",
                            "label" : "",
                            "field" : dic["internalName"],
                            "UNIPROT" : "",
                            "ENSEMBL" : dic["internalName"],
                            "add": species
                        }
                        data_list.append(data)

        with open(uniprot_schema, "r") as f:
            for line in f:
                [show, field] = line.strip().split("\t")
                if field in field_list:
                    pass
                else:
                    data = {
                        "uniprot_match" : "",
                        "displayName" : show,
                        "NCBI" : "",
                        "url" : "",
                        "selected" : "UNIPROT",
                        "uniprot_use" : "",
                        "label" : "",
                        "field" : field,
                        "UNIPROT" : show,
                        "ENSEMBL" : "",
                        "add": species
                    }
                    data_list.append(data)
        print data_list
        if len(data_list) > 0:
            self.create_db_table("sgdb_logic", data_list)

        remove_id = ['sg_gene_id', 'sg_tran_id',
                     'gene_id_1020_key', 'transcript_id_1064_key', 'translation_id_1068_key',
                     'gene_biotype', 'transcript_biotype',
                     'chromosome_name', 'start_position', 'end_position', 'strand']


        order = ['entrezgene_id', 'ensembl_gene_id', 'uniprot_gn_id', "external_gene_name"]
        ord_dic = {i:n for n,i in  enumerate(order)}
        

        with open(tran_id_file, 'r') as f:
            header = f.readline()
            tran_heads = header.strip().split("\t")

        tran_ids = [t for t in tran_heads if t not in remove_id]
        tran_ids = sorted(tran_ids,key = lambda x:ord_dic.get(x, 999))
        select_tran_ids = ['entrezgene_id', 'ensembl_gene_id', 'uniprot_gn_id', "external_gene_name"]

        with open(gene_id_file, 'r') as f:
            header = f.readline()
            gene_heads = header.strip().split("\t")

        gene_ids = [t for t in gene_heads if t not in remove_id]
        gene_ids = sorted(gene_ids,key = lambda x:ord_dic.get(x, 999))
        select_gene_ids = ['entrezgene_id', 'ensembl_gene_id', "external_gene_name"]

        genome_dic = genome_col.find_one({"genbank_assembly_accession": gb_acc})
        self.update_db_record("sgdb_genome", genome_dic["_id"],
                              species = species,
                              species_name = species_name,
                              trans_name = species_name2,
                              gene_ids = gene_ids,
                              tran_ids = tran_ids,
                              gene_sel = select_gene_ids,
                              tran_sel = select_tran_ids,
                              gene_all = gene_heads,
                              tran_all = tran_heads
        )

        self.add_gene_id(species, gene_id_file, id_index_list=select_gene_ids, text_index_list=gene_ids)
        self.add_tran_id(species, tran_id_file, id_index_list=select_tran_ids, text_index_list=tran_ids)
        

    def add_gene_id(self, species, gene_id_file, id_index_list=[], text_index_list=[]):
        """
        """
        data_list = []
        with open(gene_id_file, 'r') as f:
            for dic in csv.DictReader(f, delimiter='\t'):
                dic_clean = {k:v for k,v in dic.items() if v not in ["", "\\N"]}
                data_list.append(dic_clean)

        spe_coll = 'sgdb_{}_gene'.format(species)
        self.drop_collection(spe_coll)
        self.create_db_table(spe_coll, data_list)
        if len(id_index_list) != 0:
            id_index = [(x, ASCENDING) for x in id_index_list]
            self.create_index(spe_coll, id_index, name="all_ids")

        if len(text_index_list) != 0:
            text_index = [(x, TEXT) for x in text_index_list]
            self.create_index(spe_coll, text_index, name="all_text")

    def add_tran_id(self, species, tran_id_file, id_index_list=[], text_index_list=[]):
        """
        """

        data_list = []
        with open(tran_id_file, 'r') as f:
            for dic in csv.DictReader(f, delimiter='\t'):
                dic_clean = {k:v for k,v in dic.items() if v not in ["", "\\N"]}
                data_list.append(dic_clean)
        print len(data_list)

        spe_coll = 'sgdb_{}_tran'.format(species)
        self.drop_collection(spe_coll)
        self.create_db_table(spe_coll, data_list)
        if len(id_index_list) != 0:
            id_index = [(x, ASCENDING) for x in id_index_list]
            self.create_index(spe_coll, id_index, name="all_ids")

        if len(text_index_list) != 0:
            text_index = [(x, TEXT) for x in text_index_list]
            self.create_index(spe_coll, text_index, name="all_text")


    def search_gene_id(self, species, search_list, target_list, id_list):
        spe_coll = 'sgdb_{}_gene'.format(species)
        collect = self.db[spe_coll]
        condition = {"$or": [{k: {"$in": id_list}} for k in search_list]}
        id_result = collect.find(condition)
        return id_result


    def search_tran_id(self, species, search_list, target_list, id_list):
        spe_coll = 'sgdb_{}_tran'.format(species)
        collect = self.db[spe_coll]
        condition = {"$or": [{k: {"$in": id_list}} for k in search_list]}
        id_result = collect.find(condition)
        return id_result



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def __init__(self, method_name, params_dict=None):
        self.params_dict = params_dict
        for l in ["search_list", "target_list", "id_list"]:
            if l in self.params_dict:
                self.params_dict[l] = self.params_dict[l].split(",")
        super(TestFunction, self).__init__(methodName=method_name)
        self.toolbox = IdMapping(None)

    def add_gene_id(self):
        self.toolbox.add_gene_id(**self.params_dict)

    def add_tran_id(self):
        self.toolbox.add_tran_id(**self.params_dict)

    def add_pipe(self):
        self.toolbox.add_pipe(**self.params_dict)

    def search_gene_id(self):
        ds = self.toolbox.search_gene_id(**self.params_dict)
        for d in ds:
            print d

    def search_tran_id(self):
        ds = self.toolbox.search_tran_id(**self.params_dict)
        for d in ds:
            print d


if __name__ == '__main__':
    if sys.argv[1] in ["-h", "-help", "--h", "--help"]:
        print "\n".join(["add_gene_id", "get_tran_id", "add_pipe", "search_gene_id", "search_tran_id"])
        if len(sys.argv) == 3:
            help(getattr(IdMapping, sys.argv[2]))
            # inspect.getsourcelines(getattr(IdMapping, sys.argv[2]))
            # print h_str

    elif len(sys.argv) >= 3:
        suite = unittest.TestSuite()
        params_dict = dict()
        for par in sys.argv[2:]:
            params_dict.update({par.split("=")[0]: "=".join(par.split("=")[1:])})

        suite.addTest(TestFunction(sys.argv[1], params_dict))
        unittest.TextTestRunner(verbosity=2).run(suite)


