# !/usr/bin/python
# -*- coding: utf-8 -*-


import types
import os
import json
import unittest
import datetime
import sqlite3
import re
from collections import OrderedDict, defaultdict
from bson.son import SON
from bson.objectid import ObjectId
from mbio.api.database.prok_rna.api_base import ApiBase
from biocluster.api.database.base import Base, report_check
import glob
import copy
from pymongo import MongoClient


class ModelOrganism(ApiBase):
    def __init__(self, bind_object):
        super(ModelOrganism, self).__init__(bind_object)
        self._project_type = 'prok_rna'

    def db_update(self, database, organism=None, detail_data=None):
        update_dict = dict()
        if organism:
            update_dict['organism'] = json.dumps(organism, sort_keys=False, separators=(',', ':'))
        if detail_data:
            update_dict['detail_data'] = json.dumps(detail_data, sort_keys=False, separators=(',', ':'))
        if organism or detail_data:
            self.update_db_record('sg_model_organism_db',
                                  query_dict={'database': json.dumps(database, sort_keys=False, separators=(',', ':'))},
                                  insert_dict=update_dict)
        self.update_to_tsg_v1()

    def update_to_tsg_v1(self):
        v1_client = MongoClient("mongodb://prok_rna:y4g1v6r9d0@192.168.10.69/sanger_prok_rna?authMechanism=SCRAM-SHA-1")
        v1_db = v1_client.sanger_prok_rna
        v1_collection = v1_db["sg_model_organism_db"]
        collection = self.db["sg_model_organism_db"]
        try:
            db_dicts = collection.find()
            for each in db_dicts:
                v1_collection.insert_one(SON(each))
        except Exception:
            raise Exception("更新失败")

    def add_mo_detail(self, mo_id, result, database):
        if not os.path.exists(result):
            raise Exception('{}所指定的路径不存在，请检查！'.format(result))
        data_list = list()
        if os.path.exists(result):
            with open(result, 'r') as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    line = line.strip('\n').split('\t')
                    data = [
                        ('mo_annot_id', ObjectId(mo_id)),
                    ]
                    for i in range(len(header)):
                        if header[i] == 'coverage':
                            item = round(float(line[i])*100, 1)
                        else:
                            item = line[i]
                        data.append((header[i], item))
                    data = SON(data)
                    data_list.append(data)

            # if database == 'imodulon':
            #     detail_data = ['Gene ID', 'Gene Name', 'Gene description', 'Hit ID', 'Hit Name', 'operon', 'COG',
            #                    'Hit gene product', 'TF', 'iModulon Name', 'iM Regulator', 'iM Function', 'iM Category',
            #                    'Identity (%)', 'Coverage(%)', 'Evalue', 'Score']
            # else:
            #     detail_data = ['Gene ID', 'Gene Name', 'Gene description', 'Hit ID', 'Hit Name', 'Product',
            #                    'Transcription Factor', 'Transcription Factor Class', 'Transcription unit', 'Promoter',
            #                    'Terminator Class', 'Regulon', 'Regulon Function', 'Identity (%)', 'Coverage(%)',
            #                    'Evalue', 'Score']

            try:
                self.create_db_table('sg_model_organism_annot_detail', data_list)
                self.db['sg_model_organism_annot'].update({"_id": ObjectId(mo_id)}, {
                    "$set": {"main_id": ObjectId(mo_id), 'database': database}})
            except:
                self.bind_object.logger.error("导入模式物种注释信息：%s出错!" % (result))
            else:
                self.bind_object.logger.info("导入模式物种注释信息：%s 成功!" % result)


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
        from mbio.workflows.prok_rna.prokrna_test_api import ProkrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = dict(id='modb_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                    type='workflow',
                    name='prok_rna.prokrna_test_api',
                    options={

            },
        )
        wsheet = Sheet(data=data)
        wf = ProkrnaTestApiWorkflow(wsheet)
        wf.sheet.id = 'tsg_medical_transcriptome'
        wf.sheet.project_sn = '188_5d01dede4f911'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("prok_rna.model_organism")
        wf.test_api.db_update(
            # database={'regulon': 'RegulonDB'},
            # organism={'e_coli': 'Escherichia coli K-12'},
            # detail_data=OrderedDict([('gene_id', 'Gene ID'),
            #                          ('gene_name', 'Gene Name'),
            #                          ('description', 'Gene description'),
            #                          ('hit_id', 'Hit ID'),
            #                          ('hit_name', 'Hit Name'),
            #                          ('product', 'Product'),
            #                          ('tf', 'Transcription Factor'),
            #                          ('tf_class', 'Transcription Factor Class'),
            #                          ('tu', 'Transcription unit'),
            #                          ('promoter',  'Promoter'),
            #                          ('terminator_class', 'Terminator Class'),
            #                          ('regulon', 'Regulon'),
            #                          ('regulon_func', 'Regulon Function'),
            #                          ('identity', 'Identity (%)'),
            #                          ('coverage', 'Coverage(%)'),
            #                          ('evalue', 'Evalue'),
            #                          ('score', 'Score')]),
            database=OrderedDict([('imodulon', 'iModulonDB')]),
            organism=OrderedDict([('e_coli', 'Escherichia coli K-12'),
                                 ('s_aureus', 'Staphylococcus aureus USA 200'),
                                 ('b_subtilis', 'Bacillus subtilis')]),
            detail_data=OrderedDict([('gene_id', 'Gene ID'),
                                     ('gene_name', 'Gene Name'),
                                     ('description', 'Gene description'),
                                     ('hit_id', 'Hit ID'),
                                    ('hit_name', 'Hit Name'),
                                    ('operon','operon'),
                                    ('cog', 'COG'),
                                    ('gene_product','Hit gene product'),
                                    ('tf','TF'),
                                    ('im_name','iModulon Name'),
                                    ('im_reg','iM Regulator'),
                                    ('im_func','iM Function'),
                                    ('im_cat','iM Category'),
                                    ('identity','Identity (%)'),
                                    ('coverage','Coverage(%)'),
                                    ('evalue','Evalue'),
                                    ('score','Score')]),
        )



if __name__ == '__main__':
    unittest.main()