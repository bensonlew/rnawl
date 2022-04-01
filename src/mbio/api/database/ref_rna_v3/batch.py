# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'
from bson.objectid import ObjectId
import datetime
import glob
import json
import math
import os
import re
import types
import unittest
from collections import OrderedDict

import fastcluster as hclust
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from bson.objectid import ObjectId
from scipy import stats
from sklearn import decomposition

from mbio.api.database.whole_transcriptome.api_base import ApiBase



class Batch(ApiBase):
    def __init__(self, bind_object):
        super(Batch, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'


    def add_exp_batch(self, count_batch, exp_other, exp_id, sg_exp_batch_id,task_id,params_json):
        if type(exp_id) == str or type(exp_id) == bytes or type(exp_id) == unicode:
            exp_id = ObjectId(exp_id)
        # time_now = datetime.datetime.now()
        # name = 'exp_batch_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        #
        # # params = json.dumps({'task_id': task_id, 'submit_location': 'circrna', 'task_type': 2}, sort_keys=True)
        #
        # main_dict = {
        #     'task_id': task_id,
        #     # 'project_sn': project_sn,
        #     'name': name,
        #     'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
        #     'desc': 'exp batch main table',
        #     'params': params_json,
        #     'status': 'start',
        #     'exp_batch_id': exp_id
        # }
        # main_id = self.create_db_table('sg_exp_batch', [main_dict])
        df = pd.read_table(count_batch)
        df2 = pd.read_table(exp_other)
        exp_matrix = df.set_index('seq_id')
        other = df2.set_index('seq_id')
        exp = pd.concat([exp_matrix,other],axis=1)
        exp.index.name = 'seq_id'
        exp_final = exp.reset_index()

        exp_final['exp_id'] = exp_id
        self.create_db_table('sg_exp_detail', exp_final.to_dict('r'))
        # self.update_db_record('sg_exp_batch', sg_exp_batch_id, insert_dict={'main_id': sg_exp_batch_id})


    def add_exp_batch_corr(self, corr_output_dir, main_id=None, record_id=None):
        """
        add sample correlation result to db
        :param corr_output_dir: result of exp_corr.py
        :param exp_level:
        :param quant_method:
        :param params: which will be added to main table
        :param main_id: if provided, main table is assumed to be already created.
        :param project_sn:
        :param task_id:
        :return: main_id
        """
        results = os.listdir(corr_output_dir)
        # cluster
        samples = list()
        sample_tree = ''
        if "sample.cluster_tree.txt" in results:
            target_file = os.path.join(corr_output_dir, "sample.cluster_tree.txt")
            with open(target_file) as f:
                sample_tree = f.readline().strip()
                samples = f.readline().strip().split(";")
        else:
            print('No sample cluster result found!')
        # corr
        ordered_dict_list = list()
        if "sample_correlation.xls" in results:
            target_file = os.path.join(corr_output_dir, "sample_correlation.xls")
            corr_result = pd.read_table(target_file, index_col=0, header=0)
            corr_result = corr_result.loc[samples, :]
            # corr_result['sample'] = corr_result['sample'].astype('category')
            # corr_result['sample'].cat.reorder_categories(samples, inplace=True)
            # corr_result.sort_values('sample', inplace=True)
            corr_result.reset_index(inplace=True)
            row_dict_list = corr_result.to_dict('records')
            ordered_dict_list = self.order_row_dict_list(row_dict_list, ["sample"] + samples)
        else:
            print('No sample correlation result found!')
        # add main table
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        if not isinstance(record_id, ObjectId):
            record_id = ObjectId(record_id)
        self.update_db_record('sg_exp_batch', main_id, clust_tree=sample_tree, samples=samples, )
        tag_dict = dict(corr_id=main_id)
        self.create_db_table('sg_exp_batch_corr_detail', ordered_dict_list, tag_dict=tag_dict)
        self.update_db_record('sg_exp_batch', main_id, status="end")
        self.update_db_record('sg_exp_batch', main_id, exp_batch_id=main_id)
        self.update_db_record('sg_exp', record_id, status="end")
        return main_id, sample_tree, samples

    def add_exp_batch_pca(self, exp_output_dir, main_id, record_id):
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        if type(record_id) == str or type(record_id) == bytes or type(record_id) == unicode:
            record_id = ObjectId(record_id)
        # prepare detail table info
        # result, pc_ratio_dict = self.pca(all_exp_pd)
        target_file = os.path.join(exp_output_dir, 'PCA.xls')
        result = pd.read_table(target_file, header=0)
        row_dict_list = result.to_dict('records')
        #
        target_file = os.path.join(exp_output_dir, 'Explained_variance_ratio.xls')
        t = pd.read_table(target_file, header=None)
        pc_ratio_dict = OrderedDict(zip(t[0], t[1]))
        self.update_db_record('sg_exp_batch', main_id, ratio_dict=pc_ratio_dict)
        # insert detail
        tag_dict = dict(pca_id=main_id)
        self.create_db_table('sg_exp_batch_pca_detail', row_dict_list, tag_dict=tag_dict)
        self.update_db_record('sg_exp_batch', main_id, exp_batch_id=main_id)
        self.update_db_record('sg_exp_batch', main_id, status="end")
        return main_id

    def insert_ellipse_table(self, infile, main_id, record_id):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        if type(record_id) == str or type(record_id) == bytes or type(record_id) == unicode:
            record_id = ObjectId(record_id)
        insert_data = []
        with open(infile) as f:
            f.readline()
            name = ''
            tmp = {}
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                if spline[0] != name:
                    insert_data.append(tmp)
                    name = spline[0]
                    # name_1 = spline[0].replace('PC','').replace('_','')
                    tmp = {'name': name}
                    tmp['type'] = 'circle'
                    tmp['exp_pca_ellipse_id'] = main_id
                k = spline[1]
                number = {}
                number[k] = dict(
                    m1=spline[2],
                    c11=spline[3],
                    c12=spline[4],
                    m2=spline[5],
                    c21=spline[6],
                    c22=spline[7],
                )
                tmp[k] = json.dumps(number[k], separators=(',', ':'))
            insert_data.append(tmp)
        try:
            collection = self.db['sg_exp_batch_pca_circ_detail']
            collection.insert_many(insert_data[1:])
        except Exception as e:
            self.bind_object.set_error("导入表格sg_exp_batch_pca_circ_detail信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入sg_exp_batch_pca_circ_detail表格成功")