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
import csv
from bson.objectid import ObjectId
from api_base import ApiBase
import sys


class GenomeDb(ApiBase):
    def __init__(self, bind_object):
        super(GenomeDb, self).__init__(bind_object)

    def get_genome_info(self, species):
        genome_dict = self.db['sgdb_genome'].find_one({"species": species})
        return genome_dict

    def add_ncbi_genome_db(self, ncbi_file):
        data_list = []
        with open(ncbi_file, 'r') as f:
            for genome_dict in csv.DictReader(f, delimiter='\t'):
                genome_dict.update({"source": "NCBI"})
                genome_std = {re.sub(r'[^\da-zA-Z_]', "_", k).lower(): v for k,v in genome_dict.items()}
                data_list.append(genome_std)
        self.create_db_table("ncbi_genome_db", data_list)

    def search_ncbi(self, taxon_list):
        spe_coll = 'ncbi_genome_db'
        collect = self.db[spe_coll]
        condition = {"taxid": {"$in": taxon_list}}
        id_result = collect.find(condition)
        return id_result

    def search_ensembl(self, taxon_list):
        self.bind_object.logger.info("taxon is {}".format(taxon_list))
        spe_coll = 'ensembl_genome_db'
        collect = self.db[spe_coll]
        condition = {"taxon_id": {"$in": taxon_list}}
        id_result = collect.find(condition)
        return id_result

    def search_ensembl2(self, species_name):
        self.bind_object.logger.info("species is {}".format(species_name))
        spe_coll = 'ensembl_genome_db'
        collect = self.db[spe_coll]
        condition = {"scientific_name": species_name}
        id_result = collect.find(condition)
        return id_result

    def add_ncbi_genome_species(self, ncbi_file):
        pass

    def add_ensembl_txt(self, ensemble_file):
        data_list = []
        with open(ensemble_file, 'r') as f:
            for genome_dict in csv.DictReader(f, delimiter='\t'):
                genome_dict_new = {
                    'Accession': genome_dict.get("assembly_accession", "-"),
                    'Common name': genome_dict.get("#name", "-"),
                    'Ensembl Assembly': genome_dict.get("assembly", "-"),
                    'Genebuild Method': genome_dict.get("genebuild", "-"),
                    'Pre assembly': "-",
                    'Regulation database': "-",
                    'Scientific name': genome_dict.get("#name", "-"),
                    'Taxon ID': genome_dict.get("taxonomy_id", "-"),
                    'Variation database': genome_dict.get("variation", "-"),
                    'division': genome_dict.get("division", "-")
                }
                genome_dict_new.update({"source": "ENSEMBL"})
                genome_std = {re.sub(r'[^\da-zA-Z_]', "_", k).lower(): v for k,v in genome_dict_new.items()}
                data_list.append(genome_std)
        self.create_db_table("ensembl_genome_db", data_list)


    def add_ensembl_csv(self, ensemble_file):
        data_list = []
        with open(ensemble_file, 'r') as f:
            for genome_dict in csv.DictReader(f, delimiter=','):
                genome_dict.update({"source": "ENSEMBL",
                                    "division": "ENSEMBL"
                })
                genome_std = {re.sub(r'[^\da-zA-Z_]', "_", k).lower(): v for k,v in genome_dict.items()}
                data_list.append(genome_std)
        self.create_db_table("ensembl_genome_db", data_list)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def __init__(self, method_name, report_file):
        self.report_file = report_file
        super(TestFunction, self).__init__(methodName=method_name)
        
    toolbox = GenomeDb(None)

    def add_genome_info(self):
        self.toolbox.add_genome_info(self.report_file)

    def add_ensembl_csv(self):
        self.toolbox.add_ensembl_csv(self.report_file)

    def add_ensembl_txt(self):
        self.toolbox.add_ensembl_txt(self.report_file)

    def add_ncbi_genome_db(self):
        self.toolbox.add_ncbi_genome_db(self.report_file)



if __name__ == '__main__':
    suite = unittest.TestSuite()
    '''
    if sys.argv[1] in ["-h", "-help", "--h", "--help"]:
        print "\n".join(["add_genome_info", "get_genome_info"])
        if len(sys.argv) == 3:
    '''
    method = sys.argv[1]      
    report_file = sys.argv[2]
    suite.addTest(TestFunction(method, report_file))
    unittest.TextTestRunner(verbosity=2).run(suite)


