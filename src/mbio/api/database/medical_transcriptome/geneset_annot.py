# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:20200828

import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import pandas as pd
import json
import re
import gridfs
import glob
from collections import defaultdict
from mbio.api.database.medical_transcriptome.api_base import ApiBase
import math
import unittest


class GenesetAnnot(ApiBase):
    def __init__(self, bind_object):
        super(GenesetAnnot, self).__init__(bind_object)
        self._project_type = 'medical_transcriptome'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_main_table(self, collection_name, params, name):
        """
        添加主表的导表函数
        :param collection_name: 主表的collection名字
        :param params: 主表的参数
        :param name: 主表的名字
        :return:
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "status": "start",
            "name": name,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "params": json.dumps(params, sort_keys=True, separators=(',', ':'))
        }

        collection = self.db[collection_name]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_geneset_do_detail(self, geneset_do_table, geneset_do_id, geneset_names=None):

        """
        do详情表导表函数
        :param geneset_do_table:do结果表
        :param geneset_do_id:主表ID
        :return:
        """
        geneset_do_id = ObjectId(geneset_do_id)
        data_list = []
        geneset_names = geneset_names.split(",")
        do_df = pd.read_table(geneset_do_table, header=0)
        if do_df.shape[0]>0:
            do_df.fillna("", inplace=True)

            rename_dict = {
                "Term Type": "term_type",
                "DO Name": "do_name",
                "DO ID": "do_id"
            }
            do_df.rename(columns=rename_dict, inplace=True)
            for geneset_name in geneset_names:
                do_df[geneset_name + '_str'] = do_df[geneset_name + ' seqs']
                do_df[geneset_name + '_list'] = do_df[geneset_name + ' seqs'].map(lambda x: x.split(";") if x != "" else [])
                do_df[geneset_name + '_num'] = do_df[geneset_name + ' numbers'].map(int)
                do_df[geneset_name + '_percent'] = do_df[geneset_name + ' percent']
                do_df.drop(geneset_name + ' seqs', axis=1, inplace=True)
                do_df.drop(geneset_name + ' numbers', axis=1, inplace=True)
                do_df.drop(geneset_name + ' percent', axis=1, inplace=True)

            do_df["do_class_id"] = ObjectId(geneset_do_id)
            row_dict_list = do_df.to_dict('records')
            main_collection = self.db['sg_geneset_do_class']

            categories = list(set(do_df["term_type"]))
            try:
                self.create_db_table('sg_geneset_do_class_detail', row_dict_list)

                self.update_db_record('sg_geneset_do_class', query_dict={"_id": ObjectId(geneset_do_id)}, main_id=geneset_do_id, status="end", table_columns=geneset_names, categories=categories )

            except Exception as e:
                self.bind_object.set_error("导入main: %s出错!" , variables=("do class"))
            else:
                self.bind_object.logger.info("导入do class：%s成功!" % ("do class"))
        else:
            self.update_db_record('sg_geneset_do_class', query_dict={"_id": ObjectId(geneset_do_id)},
                                  has_results=False, status="end")

    @report_check
    def add_geneset_do_enrich(self, geneset_do_table, geneset_do_id):
        """
        do详情表导表函数
        :param geneset_do_table:do结果表
        :param geneset_do_id:主表ID
        :return:
        """
        geneset_do_id = ObjectId(geneset_do_id)
        data_list = []
        geneset_name = []
        do_df = pd.read_table(geneset_do_table, header=0)
        # 重命名
        rename_dict = dict(zip(list(do_df.columns), map(lambda x: x.strip("#").replace(" ", "_").replace("-", "_").lower(), list(do_df.columns))))
        do_df.rename(columns=rename_dict, inplace=True)
        do_df["seq_list"] = do_df["genes"].map(lambda g: g.split("|"))
        do_df["seq_str"] = do_df["genes"].map(lambda g: g.replace("|", ";"))

        do_df["do_enrich_id"] = ObjectId(geneset_do_id)

        do_df["enrich_factor"] =  do_df["ratio_in_study"].map(lambda x: float(x.split("/")[0])) /  \
                                  do_df["ratio_in_pop"].map(lambda x: float(x.split("/")[0]))

        pvalues_min =  min([p for p in do_df['pvalue'] if p > 0])
        log_pvalues_min = - math.log10(pvalues_min)
        do_df["neg_log10p_uncorrected"] = do_df['pvalue'].map(
            lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

        pvalues_min =  min([p for p in do_df['padjust'] if p > 0])
        log_pvalues_min = - math.log10(pvalues_min)
        do_df["neg_log10p_corrected"] = do_df['padjust'].map(
            lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

        categories = list(set(do_df["term_type"]))
        row_dict_list = do_df.to_dict('records')
        main_collection = self.db['sg_geneset_do_enrich']
        
        try:
            self.create_db_table('sg_geneset_do_enrich_detail', row_dict_list)
            self.update_db_record('sg_geneset_do_enrich', query_dict={"_id": geneset_do_id}, main_id=geneset_do_id, status="end", categories=categories)
            # main_collection.update({"_id": geneset_do_id}, {"main_id": geneset_do_id, "status": "end"})
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("do enrich"))
        else:
            self.bind_object.logger.info("导入do enrich ：%s出错!" % ("do enrich"))

    @report_check
    def add_geneset_reactome_enrich(self, main_id, geneset_names, geneset_reactome_table, svg_path):
        """
        reactome详情表导表函数
        :param geneset_reactome_table:reactome结果表
        :param geneset_reactome_id:主表ID
        :return:
        """
        main_id = ObjectId(main_id)
        geneset_names = geneset_names.split(",")
        reactome_df = pd.read_table(geneset_reactome_table, header=0)
        rename_dict = dict(zip(list(reactome_df.columns), map(lambda x: x.replace("#", "").replace(" ", "_").replace("-", "_").lower(), list(reactome_df.columns))))
        reactome_df.rename(columns=rename_dict, inplace=True)

        reactome_df["reactome_enrich_id"] = ObjectId(main_id)
        reactome_df["seq_list"] = reactome_df["genes"].map(lambda g: g.split("|"))
        reactome_df["seq_str"] = reactome_df["genes"].map(lambda g: g.replace("|", ";"))
        reactome_df["enrich_factor"] =  reactome_df["ratio_in_study"].map(lambda x: float(x.split("/")[0])) /  \
                                  reactome_df["ratio_in_pop"].map(lambda x: float(x.split("/")[0]))

        if len(reactome_df) > 0:
            # 考虑没有结果的情况
            pvalues_min =  min([p for p in reactome_df['pvalue'] if p > 0])
            log_pvalues_min = - math.log10(pvalues_min)
            reactome_df["neg_log10p_uncorrected"] = reactome_df['pvalue'].map(
                lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

            pvalues_min =  min([p for p in reactome_df['padjust'] if p > 0])
            log_pvalues_min = - math.log10(pvalues_min)
            reactome_df["neg_log10p_corrected"] = reactome_df['padjust'].map(
                lambda p: -math.log10(p) if p > 0 else log_pvalues_min)


        row_dict_list = reactome_df.to_dict('records')
        main_collection = self.db['sg_geneset_reactome_enrich']
        categories = list(set(reactome_df["category"]))

        self.add_geneset_reactome_enrich_svg(main_id, geneset_names, geneset_reactome_table, svg_path)

        try:
            self.create_db_table('sg_geneset_reactome_enrich_detail', row_dict_list)
            self.update_db_record('sg_geneset_reactome_enrich',
                                  query_dict={"_id": main_id},
                                  main_id=main_id, status="end",
                                  table_columns=geneset_names,
                                  categories=categories)
            # main_collection.update({"_id": geneset_reactome_id}, {"main_id": geneset_reactome_id, "status": "end"})
        except Exception as e:
            self.bind_object.logger.info("导入reactome enrich ：%s出错!" % e)
            self.bind_object.set_error("导入main: %s出错!" , variables=("reactome enrich"))
        else:
            self.bind_object.logger.info("导入reactome enrich ：%s出错!" % ("reactome enrich"))


    def add_geneset_reactome_enrich_svg(self, main_id, geneset_names, geneset_reactome_pathway, svg_path):
        data_list = []
        r_path_df = pd.read_table(geneset_reactome_pathway, header=0)
        r_path_df.fillna("", inplace=True)
        rename_dict = {
            'Pathway ID': 'pathway_id',
            'Description': 'description',
        }
        for path in r_path_df['Pathway ID']:
            with open(svg_path + '/{}.mark'.format(path), 'r') as f:
                for line in f:
                    cols = line.strip().split("\t")
                    insert_data = {
                        'reactome_enrich_id': main_id,
                        'pathway_id': path,
                        'sgbn_id': cols[2],
                        'fg_colors': cols[3],
                        'gene_list': cols[4].split(";"),
                        'geneset_list': cols[5].split(";"),
                        'title': cols[6]
                    }
                    data_list.append(insert_data)
        try:
            self.create_db_table('sg_geneset_reactome_enrich_svg', data_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("sg_geneset_reactome_enrich_pic"))
        else:
            self.bind_object.logger.info("导入：%s成功!" % ("sg_geneset_reactome_enrich_pic"))


    @report_check
    def add_geneset_reactome_detail(self, main_id, geneset_names, geneset_reactome_pathway, geneset_reactome_stat, svg_path, source=None):
        """
        reactome详情表导表函数
        :return:
        """
        data_list = []
        geneset_names = geneset_names.split(",")
        main_id = ObjectId(main_id)
        r_path_df = pd.read_table(geneset_reactome_pathway, header=0)


        self.add_geneset_reactome_pathway(main_id, geneset_names, geneset_reactome_pathway)
        self.add_geneset_reactome_stat(main_id, geneset_names, geneset_reactome_stat)
        self.add_geneset_reactome_svg(main_id, geneset_names, geneset_reactome_pathway, svg_path)

        categories = list(set(r_path_df["category"]))

        try:
            self.update_db_record('sg_geneset_reactome_class',
                                  query_dict={"_id": main_id},
                                  main_id=main_id,
                                  status="end",
                                  table_columns=geneset_names,
                                  categories=categories)

        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("sg_geneset_reactome_class"))
        else:
            self.bind_object.logger.info("导入%s成功!" % ("sg_geneset_reactome_class"))

    def add_geneset_reactome_pathway(self, main_id, geneset_names, geneset_reactome_pathway):
        r_path_df = pd.read_table(geneset_reactome_pathway, header=0)
        r_path_df.fillna("", inplace=True)

        rename_dict = {
            'Pathway ID': 'pathway_id',
            'Description': 'description',
        }
        r_path_df.rename(columns=rename_dict, inplace=True)
        r_path_df["reactome_class_id"] = main_id
        for geneset_name in geneset_names:
            r_path_df[geneset_name + '_str'] = r_path_df[geneset_name + '_genes']
            r_path_df[geneset_name + '_list'] = r_path_df[geneset_name + '_genes'].map(lambda x: x.split(";") if x != "" else [])
            r_path_df[geneset_name + '_num'] = r_path_df[geneset_name + '_numbers'].map(int)
            r_path_df.drop(geneset_name + '_genes', axis=1, inplace=True)
            r_path_df.drop(geneset_name + '_numbers', axis=1, inplace=True)

        row_dict_list = r_path_df.to_dict('records')
        try:
            self.create_db_table('sg_geneset_reactome_class_detail', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("sg_geneset_reactome_class_detail"))
        else:
            self.bind_object.logger.info("导入：%s成功!" % ("sg_geneset_reactome_class_detail"))


    def add_geneset_reactome_stat(self, main_id, geneset_names, geneset_reactome_stat):
        r_stat_df = pd.read_table(geneset_reactome_stat, header=0)
        r_stat_df.fillna("", inplace=True)
        rename_dict = {
            'Pathway ID': 'pathway_id',
            'Description': 'description',
            'Category Function description': 'category'
        }
        r_stat_df.rename(columns=rename_dict, inplace=True)
        r_stat_df["reactome_class_id"] = main_id
        for geneset_name in geneset_names:
            r_stat_df[geneset_name + '_str'] = r_stat_df[geneset_name + '_genes']
            r_stat_df[geneset_name + '_list'] = r_stat_df[geneset_name + '_genes'].map(lambda x: x.split(";") if x != "" else [])
            r_stat_df[geneset_name + '_num'] = r_stat_df[geneset_name + '_numbers'].map(int)
            r_stat_df.drop(geneset_name + '_genes', axis=1, inplace=True)
            r_stat_df.drop(geneset_name + '_numbers', axis=1, inplace=True)

        row_dict_list = r_stat_df.to_dict('records')
        try:
            self.create_db_table('sg_geneset_reactome_class_stat', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("sg_geneset_reactome_class_stat"))
        else:
            self.bind_object.logger.info("导入：%s成功!" % ("sg_geneset_reactome_class_stat"))

    def add_geneset_reactome_svg(self, main_id, geneset_name, geneset_reactome_pathway, svg_path):
        data_list = []
        r_path_df = pd.read_table(geneset_reactome_pathway, header=0)
        r_path_df.fillna("", inplace=True)
        rename_dict = {
            'Pathway ID': 'pathway_id',
            'Description': 'description',
        }
        for path in r_path_df['Pathway ID']:
            with open(svg_path + '/{}.mark'.format(path), 'r') as f:
                for line in f:
                    cols = line.strip().split("\t")
                    insert_data = {
                        'reactome_class_id': main_id,
                        'pathway_id': path,
                        'sgbn_id': cols[2],
                        'fg_colors': cols[3],
                        'gene_list': cols[4].split(";"),
                        'geneset_list': cols[5].split(";"),
                        'title': cols[6]
                    }
                    data_list.append(insert_data)
        try:
            self.create_db_table('sg_geneset_reactome_class_svg', data_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("sg_geneset_reactome_class_pic"))
        else:
            self.bind_object.logger.info("导入：%s成功!" % ("sg_geneset_reactome_class_pic"))



class TestFunction(unittest.TestCase):
    '''task
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annotation_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'medical_transcriptome.medical_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wheet)
        wf.sheet.id = 'medical_transcriptome' \
                      ''
        wf.sheet.project_sn = '188_5d12fb95c7db8'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('medical_transcriptome.geneset_annot')
        annot_merge_output_dir = '/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/AnnotMerge/output'
        params_dict = {
            'task_id': 'medical_transcriptome',
            'submit_location': 'geneset_do_class',
            'task_type': 2,
            'geneset_id': "2132131321313131,2324223232323223",
        }
        geneset_do_table = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/do_level2_class.xls"
        main_id = wf.test_api.add_main_table('sg_geneset_do_class', params_dict, "sg_geneset_do_class_test")
        wf.test_api.add_geneset_do_detail(geneset_do_table, main_id)

        params_dict = {
            'task_id': 'medical_transcriptome',
            'submit_location': 'geneset_do_class',
            'task_type': 2,
            'geneset_id': "2132131321313131,2324223232323223",
        }
        geneset_do_table = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/do_enrichment.xls"
        main_id = wf.test_api.add_main_table('sg_geneset_do_enrich', params_dict, "sg_geneset_do_enrich_test")
        wf.test_api.add_geneset_do_enrich(geneset_do_table, main_id)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

        
