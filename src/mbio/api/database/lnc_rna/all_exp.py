# !/usr/bin/python
# -*- coding: utf-8 -*-
import csv
import types
import re
import os
import json
import pandas as pd
import numpy as np
from scipy import stats
import math
from sklearn import decomposition
import scipy.cluster.hierarchy as sch
import fastcluster as hclust
from collections import OrderedDict
import unittest
import datetime
import glob
from bson.objectid import ObjectId
from api_base import ApiBase


class AllExp(ApiBase):
    def __init__(self, bind_object):
        super(AllExp, self).__init__(bind_object)

    def add_t2g_info(self, t2g_file, project_sn='lnc_rna', task_id='lnc_rna'):
        """
        add transcript and gene mapping file to database
        :param t2g_file: mapping file with two columns: transcript gene. No header line.
        :param project_sn: project number
        :param task_id: task id
        :return: main table id.
        """
        t2g_pd = pd.read_table(t2g_file, header=None)
        t2g_pd = t2g_pd.iloc[:, [0, 1]]
        t2g_pd.columns = ['transcript_id', 'gene_id']

        # add main table info
        name = "transcript2gene"
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            desc='transcript mapped to gene',
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            status='start',
        )
        main_id = self.create_db_table('sg_t2g', [main_info])

        # add detail info
        # t2g_pd["t2g_id"] = main_id
        t2g_pd["t2g_id"] = main_id
        row_dict_list = t2g_pd.to_dict('records')
        self.create_db_table('sg_t2g_detail', row_dict_list)
        self.update_db_record('sg_t2g', main_id, status="end", main_id=main_id)
        return main_id

    @staticmethod
    def process_exp_matrix(exp_matrix, log_base=None, group_dict=None):
        if type(exp_matrix) == str or type(exp_matrix) == bytes or isinstance(exp_matrix, unicode):
            all_exp_pd = pd.read_table(exp_matrix, index_col=0, header=0)
        else:
            print(exp_matrix, 'is assumed to be a pandas DataFrame Object')
            all_exp_pd = exp_matrix
        all_exp_pd.index.name = 'seq_id'

        if group_dict is not None:
            group_exp = list()
            for g in group_dict:
                g_exp = all_exp_pd.loc[:, group_dict[g]].mean(axis=1)
                g_exp.name = g
                group_exp.append(g_exp)
            all_exp_pd = pd.concat(group_exp, axis=1)

        if log_base:
            if log_base == math.e:
                all_exp_pd = np.log(all_exp_pd)
            elif log_base == 2:
                all_exp_pd = np.log2(all_exp_pd)
            elif log_base == 10:
                all_exp_pd = np.log10(all_exp_pd)
            else:
                raise Exception('log base of {} is not supported'.format(log_base))
        # return
        return all_exp_pd

    def process_exp_matrix1(exp_matrix, log_base=None, group_dict=None):
        if type(exp_matrix) == str or type(exp_matrix) == bytes or isinstance(exp_matrix, unicode):
            all_exp_pd = pd.read_table(exp_matrix, index_col=0, header=0)
        else:
            print(exp_matrix, 'is assumed to be a pandas DataFrame Object')
            all_exp_pd = exp_matrix
        all_exp_pd.index.name = 'seq_id'
        if group_dict is not None:
            group_exp = list()
            for g in group_dict:
                g_exp = all_exp_pd.loc[:, group_dict[g]].mean(axis=1)
                g_exp.name = g
                group_exp.append(g_exp)
            all_exp_pd = pd.concat(group_exp, axis=1)

        if log_base:
            if log_base == math.e:
                all_exp_pd = np.log(all_exp_pd)
            elif log_base == 2:
                all_exp_pd = np.log2(all_exp_pd)
            elif log_base == 10:
                all_exp_pd = np.log10(all_exp_pd )
            else:
                raise Exception('log base of {} is not supported'.format(log_base))
        # return
        return all_exp_pd

    def add_exp(self, exp_matrix, quant_method='RSEM', exp_level='T', group_dict=None, main_id=None, group_id=None,
                add_distribution=True, exp_type='tpm', project_sn='lnc_rna', task_id='lnc_rna',
                lib_type=None, params=None):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param exp_level: str, transcript or gene
        :param exp_type: str, usually is tpm or fpkm or count. default tpm
        :param group_dict: ordered dict of group info
        :param quant_method: the method to be used to quant expression.
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param group_id: 包含分组方案信息的主表id
        :param add_distribution: 是否添加表达分布信息到数据库
        :param params: parameters dict for expression quant.
        :return: main table id
        """

        all_exp_pd = self.process_exp_matrix(exp_matrix)
        # add main table info
        if main_id is None:
            name = "Exp" + '_' + exp_level + '_' + quant_method + '_' + exp_type.upper() + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                exp_level=exp_level,
                exp_type=exp_type.upper(),
                method=quant_method,
                desc='{} exp main table'.format(exp_level),
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                params=params,
                libtype=lib_type,
                # sample_order=sample_order + group_order,
                status="start"
            )
            main_id = self.create_db_table('sg_exp', [main_info])

        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        # update main table
        count_path = exp_matrix.replace('{}.matrix'.format(exp_type), 'count.matrix')
        #self.update_db_record('sg_exp', main_id, count_file=count_path, )
        record_dict = {"_id": main_id, "task_id": task_id}
        self.update_db_record('sg_exp', query_dict=record_dict, count_file=count_path, main_id=main_id)

        # -- add detail table info
        all_exp_pd['exp_id'] = main_id
        #all_exp_pd['is_new'] = [True if x.startswith(('MSTRG', 'TCON', 'XLOC')) else False for x in all_exp_pd.index]
        all_exp_pd.reset_index(level=0, inplace=True)
        row_dict_list = all_exp_pd.to_dict('records')
        id2gene_name, id2desc, id2gid = self.get_gene_name_dict(main_id, 'sg_exp')
        for each in row_dict_list:
           # if each['seq_id'] in id2gene_name:
                #each['gene_name'] = id2gene_name[each['seq_id']]
            #else:
                #each['gene_name'] = ''
            #each['description'] = id2desc[each['seq_id']]
            if id2gid:
                if each['seq_id'] in id2gene_name:
                    each['gene_id'] = id2gid[each['seq_id']]
                else:
                     each['gene_id']=""
             #if id2gid:
                 #each['gene_id'] = id2gid[each['seq_id']]
        self.create_db_table('sg_exp_detail', row_dict_list, )
        # add graph info
        if add_distribution:
            params = dict(
                task_id=task_id,
                exp_id=str(main_id),
                group_dict=group_dict,
                group_id=str(group_id),
                submit_location="expgraph",
                task_type=2,
                exp_level=exp_level,
                type='all',
                # quant_method=quant_method,
            )
            self.add_distribution1(exp_matrix, group_dict, params=params, exp_level=exp_level,
                                  quant_method=quant_method, project_sn=project_sn, task_id=task_id)
        # update status
        # self.update_db_record('sg_exp', main_id, status="end", main_id=main_id)
        record_dict = {"_id": main_id, "task_id": task_id}
        self.update_db_record('sg_exp', query_dict=record_dict, status="end", main_id=main_id)
        return main_id




    @staticmethod
    def get_density(all_exp_pd):
        """
        sampling 1000 density point for each columns of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        """
        records = list()
        target_columns = all_exp_pd.columns
        for sample in target_columns:
            exp = all_exp_pd[sample]
            exp = exp[exp != 0]
            try:
                density_func = stats.gaussian_kde(exp)
                min_exp, max_exp = exp.min(), exp.max()
                x_data = np.linspace(min_exp, max_exp, num=1000, endpoint=False)
                y_data = density_func(x_data)
                point_dict_list = pd.DataFrame({'log2exp': x_data, 'density': y_data}).to_dict('records')
                records.append(dict(sample=sample, data=point_dict_list))
            except Exception as e:
                print "该样本{} 不能生成密度图".format(sample)
        return records

    def add_density(self, exp_matrix, group_dict, graph_id):
        """
        Dump expression data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param group_dict: group info dict
        :param graph_id: main table id
        """
        # only use log2 transformed values for drawing
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 1]
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        records = self.get_density(all_exp_pd)
        self.create_db_table('sg_exp_graph_density', records, tag_dict=dict(graph_id=graph_id, type='sample'))
        # for group info
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        records = self.get_density(all_exp_pd)
        self.create_db_table('sg_exp_graph_density', records, tag_dict=dict(graph_id=graph_id, type='group'))

    def add_density1(self, exp_matrix, group_dict, graph_id):
        """
        Dump expression data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param group_dict: group info dict
        :param graph_id: main table id
        """
        # only use log2 transformed values for drawing
        all_exp_pd = self.process_exp_matrix1(exp_matrix, log_base=10)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 1]
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        records = self.get_density(all_exp_pd)
        self.create_db_table('sg_exp_graph_density', records, tag_dict=dict(graph_id=graph_id, type='sample'))
        # for group info
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        records = self.get_density(all_exp_pd)
        self.create_db_table('sg_exp_graph_density', records, tag_dict=dict(graph_id=graph_id, type='group'))

    @staticmethod
    def get_box(all_exp_pd):
        """
        get box plot info for each column of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        """
        stat_dict_list = list()
        target_columns = all_exp_pd.columns
        for each in target_columns:
            exp_pd = all_exp_pd[each]
            exp_pd = exp_pd[exp_pd != 0]
            summary = exp_pd.describe()
            summary.index = [u'count', u'mean', u'std', u'min', u'q1', u'median', u'q3', u'max']
            tmp_dict = summary.to_dict()
            lt25 = exp_pd[exp_pd <= tmp_dict['q1']].shape[0]
            lt50 = exp_pd[exp_pd <= tmp_dict['median']].shape[0]
            lt75 = exp_pd[exp_pd <= tmp_dict['q3']].shape[0]
            upper_whisker = tmp_dict['q3'] + 1.5*(tmp_dict['q3'] - tmp_dict['q1'])
            lower_whisker = tmp_dict['q1'] - 1.5*(tmp_dict['q3'] - tmp_dict['q1'])
            upper_outliers = list(exp_pd[exp_pd > upper_whisker])
            lower_outliers = list(exp_pd[exp_pd < lower_whisker])
            tmp_dict.update({
                'sample': each,
                'min-q1': lt25,
                'q1-median': lt50-lt25,
                'median-q3': lt75-lt50,
                'q3-max': exp_pd.shape[0]-lt75,
                'upper_whisker': upper_whisker,
                'lower_whisker': lower_whisker,
                'upper_outliers': upper_outliers,
                'lower_outliers': lower_outliers,
            })
            stat_dict_list.append(tmp_dict)
        return stat_dict_list

    def add_box(self, exp_matrix, group_dict, graph_id):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param group_dict: group dict
        :param graph_id: main table id
        """
        # only use log2 transformed values for drawing
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 1]
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        # add detail info
        stat_dict_list = self.get_box(all_exp_pd)
        self.create_db_table('sg_exp_graph_box', stat_dict_list, tag_dict=dict(graph_id=graph_id, type='sample'))
        # for group info
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        stat_list = self.get_box(all_exp_pd)
        self.create_db_table('sg_exp_graph_box', stat_list, tag_dict=dict(graph_id=graph_id, type='group'))

    def add_box1(self, exp_matrix, group_dict, graph_id):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param group_dict: group dict
        :param graph_id: main table id
        """
        # only use log2 transformed values for drawing
        all_exp_pd = self.process_exp_matrix1(exp_matrix, log_base=10)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 1]
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        # add detail info
        stat_dict_list = self.get_box(all_exp_pd)
        self.create_db_table('sg_exp_graph_box', stat_dict_list, tag_dict=dict(graph_id=graph_id, type='sample'))
        # for group info
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        stat_list = self.get_box(all_exp_pd)
        self.create_db_table('sg_exp_graph_box', stat_list, tag_dict=dict(graph_id=graph_id, type='group'))

    def add_volin(self, exp_matrix, group_dict, graph_id):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param group_dict: group dict
        :param graph_id: main table id
        """
        # only use log10 transformed values for drawing
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 1]
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        all_exp_pd = all_exp_pd.sample(frac=0.7)
        all_exp_pd.reset_index(level=0, inplace=True)
        row_dict_list = all_exp_pd.to_dict('records')
        tag_dict = dict(graph_id=graph_id, type='sample')
        self.create_db_table('sg_exp_graph_volin', row_dict_list, tag_dict=tag_dict)
        # for group info
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > -0.001]
        all_exp_pd = all_exp_pd.sample(frac=0.7)
        all_exp_pd.reset_index(level=0, inplace=True)
        row_dict_list = all_exp_pd.to_dict('records')
        tag_dict = dict(graph_id=graph_id, type='group')
        self.create_db_table('sg_exp_graph_volin', row_dict_list, tag_dict=tag_dict)

    def add_volin1(self, exp_matrix, group_dict, graph_id):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param group_dict: group dict
        :param graph_id: main table id
        """
        # only use log10 transformed values for drawing
        all_exp_pd = self.process_exp_matrix1(exp_matrix, log_base=10)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 1]
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        all_exp_pd = all_exp_pd.sample(frac=0.7)
        all_exp_pd.reset_index(level=0, inplace=True)
        row_dict_list = all_exp_pd.to_dict('records')
        tag_dict = dict(graph_id=graph_id, type='sample')
        self.create_db_table('sg_exp_graph_volin', row_dict_list, tag_dict=tag_dict)
        # for group info
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        all_exp_pd = all_exp_pd.sample(frac=0.7)
        all_exp_pd.reset_index(level=0, inplace=True)
        row_dict_list = all_exp_pd.to_dict('records')
        tag_dict = dict(graph_id=graph_id, type='group')
        self.create_db_table('sg_exp_graph_volin', row_dict_list, tag_dict=tag_dict)

    def add_distribution(self, exp_matrix, group_dict, params=None, exp_level='T',
                         quant_method='RSEM', project_sn='lnc_rna', main_id=None,
                         task_id='lnc_rna', ):
        """
        this is designed for interaction
        :param exp_matrix: exp matrix path
        :param group_dict: ordered dict
        :param params: dict of parameters for interaction
        :param quant_method: exp quant method
        :param exp_level: T or G
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :return:
        """
        if main_id is None:
            # create main table
            if type(params) == dict:
                # print(params)
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            name = "ExpDistribution" + '_' + exp_level + '_' + quant_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='density and box and volin plot main table',
                params=params,
                status="start"
            )
            main_id = self.create_db_table('sg_exp_graph', [main_info])
        else:
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        # add detail
        self.add_box(exp_matrix, group_dict, main_id)
        self.add_density(exp_matrix, group_dict, main_id)
        self.add_volin(exp_matrix, group_dict, main_id)
        self.update_db_record('sg_exp_graph', main_id, status="end", main_id=main_id)
        return main_id

    def add_distribution1(self, exp_matrix, group_dict, params=None, exp_level='T',
                         quant_method='RSEM', project_sn='lnc_rna', main_id=None,
                         task_id='lnc_rna', ):
        """
        this is designed for interaction
        :param exp_matrix: exp matrix path
        :param group_dict: ordered dict
        :param params: dict of parameters for interaction
        :param quant_method: exp quant method
        :param exp_level: T or G
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :return:
        """
        if main_id is None:
            # create main table
            if type(params) == dict:
                # print(params)
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            name = "ExpDistribution" + '_' + exp_level + '_' + quant_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='density and box and volin plot main table',
                params=params,
                status="start"
            )
            main_id = self.create_db_table('sg_exp_graph', [main_info])
        else:
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        # add detail
        self.add_box1(exp_matrix, group_dict, main_id)
        self.add_density1(exp_matrix, group_dict, main_id)
        self.add_volin1(exp_matrix, group_dict, main_id)
        self.update_db_record('sg_exp_graph', main_id, status="end", main_id=main_id)
        return main_id

    def add_venn_tt(self, venn_id):
        try:
            venn_id = ObjectId(venn_id)
        except:
            pass
        conn = self.db['sg_geneset_venn']
        conn.update({'_id': venn_id},{'$set': {'workflow': 'used'}})
        return ''

    def add_rsem_mapping(self, align_rate_file, project_sn='lnc_rna', task_id='lnc_rna'):
        "该函数暂时无用"
        align = pd.read_table(align_rate_file, index_col=0, header=0, sep='\t')
        name = 'alignment information'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            method="RSEM",
            desc='Detail of mapping reads back',
            table=align.to_dict('index'),
            status="start"
        )
        main_id = self.create_db_table('sg_alignment', [main_info])
        self.update_db_record('sg_alignment', main_id, status="end", main_id=main_id)
        return main_id

    @staticmethod
    def pca(all_exp_pd):
        """
        PCA analysis using sk-learn package of python
        :param all_exp_pd: DataFrame of expression matrix
        :return: DataFrame with samples as rows and component values as columns.
        """
        pca = decomposition.PCA()
        data = all_exp_pd.transpose()
        pca.fit(data)
        _ratio = list(enumerate(pca.explained_variance_ratio_, start=1))
        total_ratio, n_components = 0, 0
        for ind, each in _ratio:
            total_ratio += each
            if total_ratio >= 0.95:
                n_components = ind
                break
        if n_components <= 1:
            n_components = 2
        _ratio = _ratio[:n_components]
        pc_ratio = {'PC'+str(n): r for n, r in _ratio}
        result = pd.DataFrame(pca.transform(data), index=data.index)
        result = result.iloc[:, :n_components]
        result.index.name = 'sample'
        # result.columns = ['PC'+str(n)+'('+'{:.2f}%'.format(r*100)+')' for n, r in _ratio]
        result.columns = ['PC'+str(n) for n in range(1, result.shape[1] + 1)]
        return result, pc_ratio

    def add_exp_pca(self, exp_matrix, quant_method='RSEM', exp_id=None, exp_level='T', params=None,
                    project_sn='lnc_rna', task_id='lnc_rna', ):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param exp_level: str, transcript or gene
        :param exp_id: exp main table id
        :param project_sn: project id
        :param task_id: task id
        :param quant_method: method for exp quant.
        :param params: it is from POST, designed for judgement of whether the task is repeated.
        :return: main table id, one or more PCA graph data will be dumped to database.
        """
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10)
        # prepare main table info
        name = "ExpPCA" + '_' + exp_level + '_' + quant_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=exp_id,
            desc='PCA main table',
            params=params,
            status="start"
        )
        # prepare detail table info
        result, pc_ratio_dict = self.pca(all_exp_pd)
        result.reset_index(level=0, inplace=True)
        row_dict_list = result.to_dict('records')
        # insert to db
        main_info.update({"ratio_dict": pc_ratio_dict})
        main_id = self.create_db_table('sg_exp_pca', [main_info])
        tag_dict = dict(pca_id=main_id)
        self.create_db_table('sg_exp_pca_detail', row_dict_list, tag_dict=tag_dict)
        self.update_db_record('sg_exp_pca', main_id, status="end", main_id=main_id)

    def add_exp_pca2(self, exp_output_dir, quant_method='RSEM', exp_id=None, exp_level='T', params=None,
                    project_sn='lnc_rna', task_id='lnc_rna', main_id=None):
        if main_id is None:
            # prepare main table info
            name = "ExpPCA" + '_' + exp_level + '_' + quant_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                exp_id=exp_id,
                desc='PCA main table',
                params=params,
                status="start"
            )
            main_id = self.create_db_table('sg_exp_pca', [main_info])
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        # prepare detail table info
        # result, pc_ratio_dict = self.pca(all_exp_pd)
        target_file = os.path.join(exp_output_dir, 'PCA.xls')
        result = pd.read_table(target_file, header=0)
        row_dict_list = result.to_dict('records')
        #
        target_file = os.path.join(exp_output_dir, 'Explained_variance_ratio.xls')
        t = pd.read_table(target_file, header=None)
        pc_ratio_dict = OrderedDict(zip(t[0], t[1]))
        self.update_db_record('sg_exp_pca', main_id, ratio_dict=pc_ratio_dict)
        # insert detail
        tag_dict = dict(pca_id=main_id)
        self.create_db_table('sg_exp_pca_detail', row_dict_list, tag_dict=tag_dict)
        self.update_db_record('sg_exp_pca', main_id, status="end", main_id=main_id)
        return main_id

    def insert_ellipse_table(self, infile, main_id):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
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
                    #name_1 = spline[0].replace('PC','').replace('_','')
                    tmp = {'name': name}
                    tmp['type'] = 'circle'
                    tmp['exp_pca_ellipse_id'] = main_id
                k = spline[1]
                number={}
                number[k]=dict(
                    m1=spline[2],
                    c11=spline[3],
                    c12=spline[4],
                    m2=spline[5],
                    c21=spline[6],
                    c22=spline[7],
                )
                tmp[k]=json.dumps(number[k], separators=(',', ':'))
            insert_data.append(tmp)
        try:
            collection = self.db['sg_exp_pca_circ_detail']
            collection.insert_many(insert_data[1:])
        except Exception as e:
            self.bind_object.set_error("导入beta_diversity_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入beta_diversity_detail表格成功")
        main_collection = self.db['beta_diversity']
        main_collection.update_one({'_id': main_id},{'$set': {'has_ci': 1}})

    @staticmethod
    def corr(exp_pd, method='pearson'):
        """
        Correlation calculation with pandas.
        :param exp_pd: DataFrame of expression matrix
        :param method: Choices: [‘pearson’, ‘kendall’, ‘spearman’]
                        pearson : standard correlation coefficient
                        kendall : Kendall Tau correlation coefficient
                        spearman : Spearman rank correlation
        :return: correlation matrix
        """
        result = exp_pd.corr(method=method)
        result.index.name = 'sample'
        return result

    @staticmethod
    def get_hcluster_tree(exp_pd, transpose=True, method='average', metric='correlation'):
        """
        'fastcluster' was used, http://www.danifold.net/fastcluster.html?section=3.
        scipy.cluster.hierarchy.linkage could also do the same thing but slower.
        However, the documentation of scipy.cluster.hierarchy.linkage is pretty good.
        :param exp_pd: pandas DataFrame, expression matrix
        :param transpose: if to transpose the expression matrix
        :param method: methods for calculating the distance between the newly formed clusters.
            Choices: ['single', 'average', 'weighted', 'centroid', 'complete', 'median', 'ward']
        :param metric: The distance metric to use.
            Choices: ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation',
            'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski',
            'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
            'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', ]
        :return: a string denotes cluster tree with height information contained for drawing.
        """

        def make_tree(zclust, labels):
            """
            tree, subcluster, z
            :param zclust： hclust.linkage result
            :param labels: leaf label from DataFrame.columns
            :return: a string denotes cluster tree with distance information contained.
            """
            clust_step_detail = dict()
            sample_num = zclust.shape[0] + 1
            new_node = zclust.shape[0]
            for i, j, h, _ in zclust:
                i, j = int(i), int(j)
                new_node += 1
                if i < sample_num:
                    i = labels[i] + ':{:7f}'.format(h)
                if j < sample_num:
                    j = labels[j] + ':{:7f}'.format(h)
                if i in clust_step_detail:
                    i = clust_step_detail[i]
                if j in clust_step_detail:
                    j = clust_step_detail[j]
                clust_step_detail[new_node] = '({i},{j}):{h:7f}'.format(i=i, j=j, h=h)
            else:
                tree = clust_step_detail[new_node]
            # print(clust_step_detail)
            tree_list = labels[sch.leaves_list(zclust)]
            return '(' + tree + ')', tree_list
        if transpose:
            z = hclust.linkage(exp_pd.transpose(), method=method, metric=metric)
        else:
            z = hclust.linkage(exp_pd, method=method, metric=metric)
        return make_tree(z, exp_pd.columns)

    def add_exp_corr(self, exp_matrix, exp_level='T', quant_method='RSEM', sample_order=None,
                     corr_method='pearson', clust_method='average', clust_metric='euclidean',
                     project_sn='lnc_rna', task_id='lnc_rna', params=None):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param exp_level: str, transcript or gene
        :param quant_method: method for exp quant
        :param corr_method: correlation method.
        :param clust_metric: distance metric for clustering.
        :param clust_method: linkage method for clustering
        :param project_sn: project id
        :param task_id: task id.
        :param params: it is from POST, designed for judgement of whether the task is repeated.
        :return: None
        """
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10)
        # -- create main table
        name = "ExpCorr" + '_' + exp_level + '_' + quant_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='correlation analysis main table',
            params=params,
            status="start"
        )
        clust_tree, tree_list = self.get_hcluster_tree(
            all_exp_pd, transpose=True, method=clust_method, metric=clust_metric,
        )
        main_info.update({"clust_tree": clust_tree})
        main_info.update({"samples": list(tree_list)})
        main_id = self.create_db_table('sg_exp_corr', [main_info])
        # add detail
        corr_result = self.corr(all_exp_pd, method=corr_method)
        corr_result.reset_index(inplace=True)
        row_dict_list = corr_result.to_dict('records')
        ordered_dict_list = list()
        if sample_order is None:
            sample_order = tree_list
        for each in row_dict_list:
            tmp_dict = OrderedDict()
            tmp_dict["sample"] = each["sample"]
            for sample in sample_order:
                if sample in each:
                    tmp_dict[sample] = each[sample]
                else:
                    print("Sample named {} is not in sample_order".format(sample))
            ordered_dict_list.append(tmp_dict)
        tag_dict = dict(corr_id=main_id)
        self.create_db_table('sg_exp_corr_detail', ordered_dict_list, tag_dict=tag_dict)
        self.update_db_record('sg_exp_corr', main_id, status="end", main_id=main_id)

    def add_exp_corr2(self, corr_output_dir, exp_level='T', quant_method='RSEM', params=None, main_id=None,
                      project_sn='lnc_rna', task_id='lnc_rna'):
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
        if main_id is None:
            # -- create main table
            name = "ExpCorr" + '_' + exp_level + '_' + quant_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='correlation analysis main table',
                params=params,
                status="start",
            )
            main_id = self.create_db_table('sg_exp_corr', [main_info])
        # add detail table
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        self.update_db_record('sg_exp_corr', main_id, clust_tree=sample_tree, samples=samples,)
        tag_dict = dict(corr_id=main_id)
        self.create_db_table('sg_exp_corr_detail', ordered_dict_list, tag_dict=tag_dict)
        self.update_db_record('sg_exp_corr', main_id, status="end", main_id=main_id)
        return main_id, sample_tree, samples

    def add_exp_venn(self, venn_graph, venn_table=None, project_sn='lnc_rna', exp_level='T', main_id=None,
                     quant_method='RSEM', task_id='lnc_rna', params=None):
        """
        add venn analysis info
        :param venn_graph: venn_graph.xls resulted from express_venn tool
        :param venn_table: venn_table.xls resulted from express_venn tool
        :param quant_method: exp quant method
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param exp_level: T or G
        :param params: string of parameters dict, designed for judgement of whether the task is repeated.
        :return: main table id
        """
        # add main table info
        if main_id is None:
            name = "ExpVenn" + '_' + exp_level + '_' + quant_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='venn main table',
                params=params,
                status="start"
            )
            main_id = self.create_db_table('sg_exp_venn', [main_info])
        else:
            main_id = ObjectId(main_id)
        # add detail table info
        graph_pd = pd.read_table(venn_graph, header=0, sep='\t', keep_default_na = False)
        graph_pd.columns = ["name", "seqs"]
        detail_dict_list = graph_pd.to_dict('records')
        if venn_table:
            table_pd = pd.read_table(venn_table, header=None, sep='\t', keep_default_na = False)
            table_pd.columns = ["combination", "num", "only_list"]
            detail_dict_list += table_pd.to_dict('records')
        self.create_db_table('sg_exp_venn_detail', detail_dict_list, tag_dict={'venn_id': main_id})
        self.update_db_record('sg_exp_venn', main_id, status="end", main_id=main_id)
        return main_id

    @staticmethod
    def _get_volcano_status_cutoff(diff_table, pvalue_padjust):
        sig_status = list()
        sig_mark = diff_table['significant']
        reg_list = diff_table['regulate']
        if 'no' in list(sig_mark):
            no_sig_num = sig_mark[sig_mark == "no"].shape[0]
            sig_status.append('nosig_'+str(no_sig_num))
        if 'yes' in list(sig_mark):
            reg_mark = reg_list[sig_mark == 'yes']
            if 'down' in list(reg_mark):
                down_num = reg_mark[reg_mark == 'down'].shape[0]
                sig_status.append('down_'+str(down_num))
            if 'up' in list(reg_mark):
                up_num = reg_mark[reg_mark == 'up'].shape[0]
                sig_status.append('up_'+str(up_num))

        sig_pvalues = diff_table[pvalue_padjust][diff_table['significant'] == "yes"]
        log10_sig_pvalues = -np.log10(sig_pvalues)
        log10_pvalue_list = sorted(list(log10_sig_pvalues[log10_sig_pvalues > 0]))

        if len(sig_pvalues) > 2000:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.85)]
        elif len(sig_pvalues) > 1000:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.90)]
        elif len(sig_pvalues) > 500:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.95)]
        elif len(sig_pvalues) > 250:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.99)]
        elif len(sig_pvalues) == 0:
            tmp = -np.log10(diff_table[pvalue_padjust])
            tmp_list = sorted(tmp[tmp > 0])
            if len(tmp_list) == 0:
                log10_pvalue_cutoff = 200
            else:
                log10_pvalue_cutoff = tmp_list[int(len(tmp_list)*0.9)]
        else:
            # print(pvalue_padjust, diff_table, log10_pvalue_list)
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.8)]
        return sig_status, log10_pvalue_cutoff

    def add_diffexp(self, diff_output, exp_id=None, group_dict=None, group_id=None,
                    exp_level='T', quant_method='RSEM', diff_method='DESeq2', main_id=None,
                    project_sn='lnc_rna', task_id='lnc_rna', params=None,rna_type="all",
                    pvalue_padjust='padjust', create_geneset=True, seq_type=None):
        """
        add differential analysis result to database
        :param diff_output: diffexp result dir
        :param exp_id: exp table id from POST, will be used for getting expression value
        :param group_dict: group info dict. {group:[s1,s2,], ...}
        :param exp_level: expression level
        :param pvalue_padjust: pvalue or padjust, for significant judgement.
        :param quant_method: method for expression quant
        :param diff_method: differential analysis method
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param group_id: 包含分组方案信息的主表id
        :param create_geneset: 是否创建差异基因集
        :param params: from POST, designed for judgement of whether the task is repeated.
        :return: main table id
        """
        summary = glob.glob(os.path.join(diff_output, '*_diff_summary.xls'))[0]
        # diff_method = os.path.basename(summary).split('_diff_summary.xls')[0]
        diff_files = glob.glob(os.path.join(diff_output, '*_vs_*.*.xls'))
        if not diff_files:
            raise Exception('No target file found in {}'.format(diff_output))
        cmp_list = list()
        cmp_detail_dict = dict()
        diff_dict_list = list()
        volcano_dict_list = list()
        sig_status = dict()
        scatter_dict_list = list()
        if seq_type and os.path.exists(seq_type):
            with open(seq_type, 'r') as f:
                mrna_list = [line.split("\t")[0] for line in f.readlines() if line.split("\t")[2] == "mRNA"]
            with open(seq_type, 'r') as f:
                lncrna_list = [line.split("\t")[0] for line in f.readlines() if line.split("\t")[2] == "lncRNA"]

        for each in diff_files:
            if each.endswith('.annot.xls') or each.endswith('.normalize.xls') or each.endswith('.sizeFactor.xls') or each.endswith('.normFactor.xls'):
                continue
            diff_pd = pd.read_table(each, header=0, sep='\t')
            columns = diff_pd.columns
            log2fc_ind = list(columns).index('log2fc')
            need_cols = ['seq_id', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
            # 分组各组表达均值
            need_cols += [columns[log2fc_ind-3], columns[log2fc_ind-2]]
            print(need_cols)
            samples = list()
            for x in columns:
                _m = re.match(r'(.*)_count$', x)
                if _m:
                    samples.append(_m.groups()[0])
            fname = os.path.basename(each)
            ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(diff_method.lower()), fname).groups()
            cmp_combine = ctrl + '|' + test
            cmp_list.append(cmp_combine)
            cmp_detail_dict[cmp_combine] = samples
            cmp_pd = pd.DataFrame([cmp_combine]*diff_pd.shape[0], columns=['compare'])
            tmp_pd = pd.concat([diff_pd.loc[:, need_cols], cmp_pd], axis=1)
            tmp_pd.columns = list(tmp_pd.columns[:-3]) + ['group1', 'group2', 'compare']
            print(tmp_pd.columns)
            diff_dict_list += tmp_pd.to_dict('records')
            # get volcano data
            status_list, stat_cutoff = self._get_volcano_status_cutoff(diff_pd, pvalue_padjust)
            sig_status[cmp_combine] = status_list
            volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', pvalue_padjust, 'significant', 'regulate']]
            bool_ind = volcano_pd[pvalue_padjust] <= 0
            min_pvalue = min([x if x>0 else '' for x in volcano_pd[pvalue_padjust].tolist()])
            volcano_pd.loc[bool_ind, pvalue_padjust] = min_pvalue
            volcano_pd[pvalue_padjust] = -np.log10(volcano_pd[pvalue_padjust])
            volcano_pd.dropna(inplace=True)
            volcano_pd.columns = ['seq_id', 'log2fc', 'log10pvalue', 'significant', 'regulate']
            bool_ind = volcano_pd['log10pvalue'] > stat_cutoff
            volcano_pd.loc[bool_ind, 'log10pvalue'] = stat_cutoff
            volcano_pd = pd.concat([volcano_pd, cmp_pd], axis=1)
            volcano_pd_nosig = volcano_pd[volcano_pd['significant'] == 'no']
            # random select 10000 not sig diff genes for plotting
            if volcano_pd_nosig.shape[0] > 8000:
                volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
                volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.5)
                if volcano_pd_nosig.shape[0] > 12000:
                    volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
                if volcano_pd_nosig.shape[0] > 12000:
                    volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
                volcano_pd = pd.concat([volcano_pd_sig, volcano_pd_nosig], axis=0)
            volcano_dict_list += volcano_pd.to_dict('records')
            # get scatter data
            scatter_pd = tmp_pd.loc[:, ['seq_id', 'group1', 'group2', 'compare', 'significant', 'regulate']]
            scatter_pd.set_index('seq_id', inplace=True)
            scatter_pd = scatter_pd.loc[volcano_pd['seq_id'], :].reset_index()
            scatter_pd['group1'] = (scatter_pd['group1']+1).apply(np.log10)
            scatter_pd['group2'] = (scatter_pd['group2']+1).apply(np.log10)
            scatter_dict_list += scatter_pd.to_dict('records')
            # get significant diff gene set and add to database
            if create_geneset:
                
                name = ctrl + '_vs_' + test + '_' + exp_level

                geneset_main_info = dict(
                    project_sn=project_sn,
                    task_id=task_id,
                    name=name,
                    type=exp_level,
                    desc='differential expressed gene set',
                    group_id=group_id,
                    source="diff_exp",
                    gene_length=0,
                    is_use=1
                )
                mrna_main_info = geneset_main_info.copy()
                lncrna_main_info = geneset_main_info.copy()

                if seq_type and os.path.exists(seq_type):
                    pd_select = diff_pd[(diff_pd["significant"] == "yes") & (diff_pd["seq_id"].isin(mrna_list))]
                    pd_m_list = list(pd_select["seq_id"])
                    pd_mr_list = list(pd_select["regulate"])
                    if len(pd_m_list) >= 1:
                        mrna_main_info["gene_length"] = len(pd_m_list)
                        genet_detail_info = [{"seq_list": pd_m_list, "regulate_list": pd_mr_list}]
                        self.add_set(mrna_main_info, genet_detail_info)
                    pd_select = diff_pd[(diff_pd["significant"] == "yes") & (diff_pd["seq_id"].isin(lncrna_list))]
                    pd_l_list = list(pd_select["seq_id"])
                    pd_lr_list = list(pd_select["regulate"])

                    lncrna_main_info["type"] = "L" + geneset_main_info["type"]
                    if len(pd_l_list) >= 1:
                        lncrna_main_info["gene_length"] = len(pd_l_list)
                        lncrna_main_info["name"] = ctrl + '_vs_' + test + '_L' + exp_level
                        genet_detail_info = [{"seq_list": pd_l_list, "regulate_list": pd_lr_list}]
                        self.add_set(lncrna_main_info, genet_detail_info)
                else:
                    sig_seqs = list(diff_pd['seq_id'][diff_pd['significant'] == 'yes'])
                    sig_seqs_r = list(diff_pd['regulate'][diff_pd['significant'] == 'yes'])
                    geneset_main_info["gene_length"] = len(sig_seqs)

                    if len(sig_seqs) >= 1:
                        genet_detail_info = [{"seq_list": sig_seqs, "regulate_list": sig_seqs_r}]
                        self.add_set(geneset_main_info, genet_detail_info)
        else:
            'loop end of diff_files'

        if main_id is None:
            # add main table info
            name = "Diff" + '_' + exp_level + '_' + quant_method + '_' + diff_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='differential analysis main table',
                exp_level=exp_level,
                group_dict=group_dict,
                group_id=group_id,
                rna_type=rna_type,
                # control_id=cmp_id,
                exp_id=exp_id,
                # diff_sum=dict(zip(levels[0][labels[0]], levels[1][labels[1]])),
                params=params,
                status="start",
            )
            main_id = self.create_db_table('sg_diff', [main_info])
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        # update main table
        self.update_db_record('sg_diff', main_id,
                              cmp_list=cmp_list,
                              cmp_detail=cmp_detail_dict,
                              sig_status=sig_status,
                              main_id=main_id
                              )
        # get id2gene_name dict
        id2gene_name, id2desc, id2gid = self.get_gene_name_dict(main_id, 'sg_diff')

        # add detail info
        for each in diff_dict_list:
            #each['gene_name'] = id2gene_name[each['seq_id']]
            #each['description'] = id2desc[each['seq_id']]
            if id2gid:
                each['gene_id'] = id2gid[each['seq_id']]
            if each['seq_id'] in mrna_list:
                each['rna_type']="mRNA"
            else:
                each['rna_type']="lncRNA"
        self.create_db_table('sg_diff_detail', diff_dict_list, tag_dict={'diff_id': main_id})

        # add summary detail
        with open(summary, "r") as f:
            lines = f.readlines()
            if len(lines) >= 3:
                summary_pd = pd.read_table(summary, header=[0, 1])
                levels = summary_pd.columns.levels
                labels = summary_pd.columns.labels
                summary_pd.columns = levels[0][labels[0]]
                #summary_pd['gene_name'] = [id2gene_name[x] for x in summary_pd['seq_id']]
                #summary_pd['description'] = [id2desc[x] for x in summary_pd['seq_id']]
#                if id2gid:
  #                  summary_pd['gene_id'] = [id2gid[x] for x in summary_pd['seq_id']]
                summary_dict_list = summary_pd.to_dict('records')
                for each in summary_dict_list:
                    if each['seq_id'] in mrna_list:
                        each['rna_type'] = "mRNA"
                    else:
                        each['rna_type'] = "lncRNA"
                self.create_db_table('sg_diff_summary', summary_dict_list, tag_dict={'diff_id': main_id})

                # add volcano detail
 #               for each in volcano_dict_list:
  #                  each['gene_name'] = id2gene_name[each['seq_id']]
                for each in volcano_dict_list:
                    if each['seq_id'] in mrna_list:
                        each['rna_type'] = "mRNA"
                    else:
                        each['rna_type'] = "lncRNA"
                self.create_db_table('sg_diff_volcano', volcano_dict_list, tag_dict={'diff_id': main_id})

                # add scatter detail
   #             for each in scatter_dict_list:
    #                each['gene_name'] = id2gene_name[each['seq_id']]
                for each in scatter_dict_list:
                    if each['seq_id'] in mrna_list:
                        each['rna_type'] = "mRNA"
                    else:
                        each['rna_type'] = "lncRNA"
                self.create_db_table('sg_diff_scatter', scatter_dict_list, tag_dict={'diff_id': main_id})
            else:
                pass

        # update status
        self.update_db_record('sg_diff', main_id, status="end", main_id=main_id, )
        return main_id

    def add_diffexp_all(self, uniform_output, diff_output, exp_id=None, group_dict=None, group_id=None,
                        exp_level='T', quant_method='RSEM', diff_method='DESeq2', main_id=None,
                        project_sn='lnc_rna', task_id='lnc_rna', params=None, rna_type="all",
                        pvalue_padjust='padjust', create_geneset=True, nosig_info=None, seq_type=None):
        """
        add differential analysis result to database
        :param diff_output: diffexp result dir
        :param exp_id: exp table id from POST, will be used for getting expression value
        :param group_dict: group info dict. {group:[s1,s2,], ...}
        :param exp_level: expression level
        :param pvalue_padjust: pvalue or padjust, for significant judgement.
        :param quant_method: method for expression quant
        :param diff_method: differential analysis method
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param group_id: 包含分组方案信息的主表id
        :param create_geneset: 是否创建差异基因集
        :param params: from POST, designed for judgement of whether the task is repeated.
        :return: main table id
        """
        if nosig_info:
            nosig_list = dict()
            nosig_pd = pd.read_table(nosig_info, header=0, sep='\t')
            nosig_list = nosig_pd.to_dict('records')
            self.create_db_table('sg_diff_detail', nosig_list, tag_dict={'diff_id': main_id, "seq_type": "no_sig"})

        if seq_type and os.path.exists(seq_type):
            with open(seq_type, 'r') as f:
                mrna_list = [line.split("\t")[0] for line in f.readlines() if line.split("\t")[2] == "mRNA"]
            with open(seq_type, 'r') as f:
                lncrna_list = [line.split("\t")[0] for line in f.readlines() if line.split("\t")[2] == "lncRNA"]

        summary = glob.glob(os.path.join(diff_output, '*diff_summary*.xls'))[0]
        json_path = os.path.join(uniform_output, 'json')
        js = open(json_path, 'r+')
        js_load = json.load(js)
        cmp_list = js_load['cmp_list']
        cmp_detail_dict = js_load['cmp_detail_dict']
        sig_status = js_load['sig_status']
        if main_id is None:
            # add main table info
            name = "DiffExpress" + '_' + exp_level + '_' + quant_method + '_' + diff_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1.1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='differential analysis main table',
                exp_level=exp_level,
                group_dict=group_dict,
                group_id=group_id,
                exp_id=exp_id,
                rna_type=rna_type,
                params=params,
                status="start",
            )
            main_id = self.create_db_table('sg_diff', [main_info])
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        # update main table
        self.update_db_record('sg_diff', main_id,
                              cmp_list=cmp_list,
                              cmp_detail=cmp_detail_dict,
                              sig_status=sig_status,
                              main_id=main_id,
                              has_regulate="yes"
                              )
        all_detail = pd.read_table(os.path.join(uniform_output, 'all_detail.txt'), header=0, sep='\t')
        all_detail.fillna('', inplace=True)
        diff_dict_list = all_detail.to_dict('records')
        id2gene_name, id2desc, id2gid = self.get_gene_name_dict(main_id, 'sg_diff')
        for each in diff_dict_list:
            if id2gid:
                each['gene_id'] = id2gid[each['seq_id']]
            if each['seq_id'] in mrna_list:
                each['rna_type'] = "mRNA"
            else:
                each['rna_type'] = "lncRNA"
        try:
            self.create_db_table('sg_diff_detail', diff_dict_list, tag_dict={'diff_id': main_id})
        except Exception as e:
            self.bind_object.set_error("导入表格sg_diff_detail信息出错:{}".format(e))
        if create_geneset:
            for i in all_detail.groupby('compare'):
                compare = i[0]
                diff_pd = i[1]
                ctrl, test = compare.split('|')
                name = ctrl + '_vs_' + test + '_' + exp_level
                geneset_main_info = dict(
                    project_sn=project_sn,
                    task_id=task_id,
                    name=name,
                    type=exp_level,
                    desc='differential expressed gene set',
                    group_id=group_id,
                    source="diff_exp",
                    gene_length=0,
                    is_use=1
                )
                mrna_main_info = geneset_main_info.copy()
                lncrna_main_info = geneset_main_info.copy()

                if seq_type and os.path.exists(seq_type):
                    pd_select = diff_pd[(diff_pd["significant"] == "yes") & (diff_pd["seq_id"].isin(mrna_list))]
                    pd_m_list = list(pd_select["seq_id"])
                    pd_mr_list = list(pd_select["regulate"])
                    if len(pd_m_list) >= 1:
                        mrna_main_info["gene_length"] = len(pd_m_list)
                        genet_detail_info = [{"seq_list": pd_m_list, "regulate_list": pd_mr_list}]
                        self.add_set(mrna_main_info, genet_detail_info)
                    pd_select = diff_pd[(diff_pd["significant"] == "yes") & (diff_pd["seq_id"].isin(lncrna_list))]
                    pd_l_list = list(pd_select["seq_id"])
                    pd_lr_list = list(pd_select["regulate"])

                    lncrna_main_info["type"] = "L" + geneset_main_info["type"]
                    if len(pd_l_list) >= 1:
                        lncrna_main_info["gene_length"] = len(pd_l_list)
                        lncrna_main_info["name"] = ctrl + '_vs_' + test + '_L' + exp_level
                        genet_detail_info = [{"seq_list": pd_l_list, "regulate_list": pd_lr_list}]
                        self.add_set(lncrna_main_info, genet_detail_info)
                else:
                    sig_seqs = list(diff_pd['seq_id'][diff_pd['significant'] == 'yes'])
                    sig_seqs_r = list(diff_pd['regulate'][diff_pd['significant'] == 'yes'])
                    geneset_main_info["gene_length"] = len(sig_seqs)

                    if len(sig_seqs) >= 1:
                        genet_detail_info = [{"seq_list": sig_seqs, "regulate_list": sig_seqs_r}]
                        self.add_set(geneset_main_info, genet_detail_info)
        del all_detail
        del diff_dict_list
        # add summary detail
        with open(summary, "r") as f:
            lines = f.readlines()
            if len(lines) >= 3:
                summary_pd = pd.read_table(summary, header=[0, 1])
                levels = summary_pd.columns.levels
                labels = summary_pd.columns.labels
                summary_pd.columns = levels[0][labels[0]]
                # add by fwy 20200729 将yes/no和up/down分开导入
                columns = summary_pd.columns
                for column in columns[1:-1]:
                    summary_pd[column + "_regulate"] = summary_pd.apply(lambda x: x[column].split("|")[-1], axis=1)
                    summary_pd[column] = summary_pd.apply(lambda x: x[column].split("|")[0], axis=1)

                summary_dict_list = summary_pd.to_dict('records')
                for each in summary_dict_list:
                    if each['seq_id'] in mrna_list:
                        each['rna_type'] = "mRNA"
                    else:
                        each['rna_type'] = "lncRNA"
                try:
                    self.create_db_table('sg_diff_summary', summary_dict_list, tag_dict={'diff_id': main_id})
                except Exception as e:
                    self.bind_object.set_error("导入表格sg_diff_summary信息出错:{}".format(e))
                del summary_pd
                del summary_dict_list
                all_volcano = pd.read_table(os.path.join(uniform_output, 'all_volcano.txt'), header=0, sep='\t')
                volcano_dict_list = all_volcano.to_dict('records')
                for each in volcano_dict_list:
                    if each['seq_id'] in mrna_list:
                        each['rna_type'] = "mRNA"
                    else:
                        each['rna_type'] = "lncRNA"
                try:
                    self.create_db_table('sg_diff_volcano', volcano_dict_list, tag_dict={'diff_id': main_id})
                except Exception as e:
                    self.bind_object.set_error("导入表格sg_diff_volcano信息出错:{}".format(e))
                del all_volcano
                del volcano_dict_list
                all_scatter = pd.read_table(os.path.join(uniform_output, 'all_scatter.txt'), header=0, sep='\t')
                scatter_dict_list = all_scatter.to_dict('records')
                for each in scatter_dict_list:
                    if each['seq_id'] in mrna_list:
                        each['rna_type'] = "mRNA"
                    else:
                        each['rna_type'] = "lncRNA"
                try:
                    self.create_db_table('sg_diff_scatter', scatter_dict_list, tag_dict={'diff_id': main_id})
                except Exception as e:
                    self.bind_object.set_error("导入表格sg_diff_scatter信息出错:{}".format(e))
                del all_scatter
                del scatter_dict_list
            else:
                pass

        # update status
        self.update_db_record('sg_diff', main_id, status="end", main_id=main_id, )
        return main_id

    def add_set(self, main_info, detail_info):
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        main_info.update(dict(status="start", created_ts=created_ts))
        main_id = self.create_db_table('sg_geneset', [main_info])
        self.create_db_table('sg_geneset_detail', detail_info, tag_dict={"geneset_id": main_id})
        # self.update_db_record('sg_geneset', main_id, status="end", is_use=0, main_id=main_id)
        task_id = main_info['task_id']
        record_dict = {"_id": main_id, "task_id": task_id}
        self.update_db_record('sg_geneset', query_dict=record_dict, status="end", is_use=0, main_id=main_id, params=task_id)
        return main_id

    def add_geneset_cluster(self, cluster_output_dir, gene_detail_file, main_id=None, project_sn='lnc_rna',
                            task_id='lnc_rna', params=None):
        # prepare main_table data
        results = os.listdir(cluster_output_dir)
        gene_cluster, sample_cluster = False, False
        genes, samples = list(), list()
        gene_tree, sample_tree = "", ""

        if "seq.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "seq.cluster_tree.txt")
            with open(target_file) as f:
                gene_cluster = True
                gene_tree = f.readline().strip()
                genes = f.readline().strip().split(";")
        #
        if "sample.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "sample.cluster_tree.txt")
            with open(target_file) as f:
                sample_cluster = True
                sample_tree = f.readline().strip()
                samples = f.readline().strip().split(";")
        #
        if "seq.kmeans_cluster.txt" in results:
            gene_cluster = True
            target_file = os.path.join(cluster_output_dir, "seq.kmeans_cluster.txt")
            with open(target_file) as f:
                genes = list()
                for line in f:
                    if not line.strip():
                        continue
                    genes += line.strip().split('\t')[1].split(";")
        #
        detail_info = list()
        trend_dict = dict()
        seq_id2gene_info = dict()
        with open(gene_detail_file, 'rb') as f:
            # zhaozhipeng -- 20190319 修改
            for dic in csv.DictReader(f, delimiter='\t'):
                transcript_id = dic.pop('transcript_id')
                seq_id2gene_info[transcript_id] = dic
            # for linen in f.readlines()[1:]:
            # line = linen.strip()
            # if len(line.split("\t")) > 2:
            #     if line.split("\t")[2].strip() in ["-", "_"]:
            #         seq_id2name[line.split("\t")[0]] = line.split("\t")[0]
            #         seq_id2name[line.split("\t")[1]] = line.split("\t")[1]
            #     else:
            #         seq_id2name[line.split("\t")[0]] = line.split("\t")[2].strip()
            #         seq_id2name[line.split("\t")[1]] = line.split("\t")[2].strip()
            # else:
            #     seq_id2name[line.split("\t")[0]] = line.split("\t")[0]
            #     seq_id2name[line.split("\t")[1]] = line.split("\t")[1]

        def add_col2df(m_df):
            """ zhaozhipeng -- 20190319 添加

            :param m_df:
            :return:
            """
            tmp_dic = {}
            tmp_list = ('gene_name', 'gene_id', 'gene_desc')
            gene_names, gene_ids, gene_descs = zip(*[
                [dic.get(i, '') for i in tmp_list]
                for dic in (seq_id2gene_info.get(x, tmp_dic) for x in m_df['seq_id'])
            ])
            m_df['gene_name'] = gene_names
            m_df['gene_id'] = gene_ids
            m_df['gene_desc'] = gene_descs
            # m_df['gene_desc'] = m_df['seq_id'].map(
            #     lambda x: seq_id2gene_info.get(x, tmp_dic).get('gene_desc', ''))

        if ("seq.cluster_tree.txt" in results) or ("seq.kmeans_cluster.txt" in results):
            sub_clusters = [x for x in results if x.startswith('seq.subcluster')]
            number_order = [(x, int(x.split('_')[1])) for x in sub_clusters]
            tmp = sorted(number_order, key=lambda x: x[1])
            sub_clusters = [x[0] for x in tmp]

            for sub in sub_clusters:
                target_file = os.path.join(cluster_output_dir, sub)
                tmp_df = pd.read_table(target_file, header=0)
                sub_cluster_id = int(sub.split('_')[1])
                tmp_df["sub_cluster"] = sub_cluster_id
                # zhaozhipeng -- 20190319 修改
                add_col2df(tmp_df)
                detail_info += json.loads(tmp_df.to_json(orient="records"))
                mean_dict = tmp_df.iloc[:, 1:-1].mean().to_dict()
                trend_dict[str(sub_cluster_id)] = mean_dict
        #
        target_file = os.path.join(cluster_output_dir, "expression_matrix.xls")

        exp_pd = pd.read_table(target_file, header=0)

        if not detail_info:
            # zhaozhipeng -- 20190319 修改
            add_col2df(exp_pd)
            detail_info = exp_pd.to_dict('records')

        if not genes:
            genes = list(exp_pd['seq_id'])
        if not samples:
            samples = list(exp_pd.columns)[1:]
        # add main table info'
        if main_id is None:
            name = "GeneSet_Cluster" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if params is None:
                params_dict = dict()
            elif type(params) == dict:
                params_dict = params
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            else:
                params_dict = json.loads(params)

            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='geneset cluster main table',
                status="start",
                params=params,
                type='T' if "exp_level" not in params_dict else params_dict["exp_level"],
            )
            main_id = self.create_db_table('sg_geneset_cluster', [main_info])
        else:
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
        # update main table
        self.update_db_record('sg_geneset_cluster', main_id,
                              trend_dict=trend_dict,
                              samples=samples,
                              gene_cluster=gene_cluster,
                              sample_cluster=sample_cluster, )
        # add detail info
        tree_info = dict(
            genes=genes,
            gene_tree=gene_tree,
            sample_tree=sample_tree,
            cluster_id=main_id,
        )
        self.create_db_table('sg_geneset_cluster_tree', [tree_info])
        self.create_db_table('sg_geneset_cluster_detail', detail_info, tag_dict=dict(cluster_id=main_id))
        self.update_db_record('sg_geneset_cluster', main_id, status="end", main_id=main_id, )
        return main_id

    def get_gene_name_dict(self, main_id, main_table_name):
        """为了给差异分析结果添加gene_name"""
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        diff_main = self.db[main_table_name].find_one({"main_id": main_id})
        task_id = diff_main['task_id']
        exp_level = diff_main['exp_level']
        annot_table = self.db['sg_annotation_query']
        annot_main = annot_table.find_one({"task_id": task_id, "type": "latest"})
        if not annot_main:
            annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
        if "main_id" not in annot_main:
            annot_main_id = annot_main['_id']
        else:
            annot_main_id = annot_main['main_id']
        annot_detail = self.db['sg_annotation_query_detail']
        query_dict = dict(query_id=annot_main_id, )
        result_dict = dict(_id=0, gene_name=1, gene_id=1, transcript_id=1, description=1)
        result = annot_detail.find(query_dict, result_dict)
        gene_annot = pd.DataFrame(list(result))
        if exp_level[0].lower() == 't':
            id2gene_name = dict(zip(gene_annot['transcript_id'], [x if x else '-' for x in gene_annot['gene_name']]))
            id2desc = dict(zip(gene_annot['transcript_id'], [x if x else '-' for x in gene_annot['description']]))
            id2gid = dict(zip(gene_annot['transcript_id'], [x if x else '-' for x in gene_annot['gene_id']]))
        else:
            id2gene_name = dict(zip(gene_annot['gene_id'], [x if x else '-' for x in gene_annot['gene_name']]))
            id2desc = dict(zip(gene_annot['gene_id'], [x if x else '-' for x in gene_annot['description']]))
            id2gid = dict()
        return id2gene_name, id2desc, id2gid

    def run1(self):
        self.add_exp("/mnt/ilustre/users/sanger-dev/workspace/20190415/Single_classify_quant12-41-32/ClassifyQuant/output/transcript.tpm.matrix", quant_method="salmon", exp_level='T', group_dict=None, main_id=None, group_id=None,add_distribution=False, exp_type='tpm', project_sn='lnc_rna', task_id='lnc_rna',lib_type=None, params=None)
        self.add_exp("/mnt/ilustre/users/sanger-dev/workspace/20190415/Single_classify_quant12-41-32/ClassifyQuant/output/gene.tpm.matrix",quant_method="salmon", exp_level='G', group_dict=None, main_id=None, group_id=None,add_distribution=False, exp_type='tpm', project_sn='lnc_rna', task_id='lnc_rna',lib_type=None, params=None)

    def run2(self):
        self.add_t2g_info("/mnt/ilustre/users/sanger-dev/workspace/20190322/Single_assemble_3643_9041/Assemble/output/NewTranscripts/trans2gene")
class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_mongo(test):
      from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
      from biocluster.wsheet import Sheet
      import random

      data = {
        # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
        "id": "lnc_rna_"+str(random.randint(1,10000)),
        "project_sn": "lnc_rna",
        "type": "workflow",
        "name": "lnc_rna.lnc_rna_test_api",
        "options": {
        },
      }
      wsheet = Sheet(data=data)
      wf = LncRnaTestApiWorkflow(wsheet)
      wf.IMPORT_REPORT_DATA = True
      wf.IMPORT_REPORT_AFTER_END = False
      wf.test_api = wf.api.api("lnc_rna.all_exp")
      wf.test_api.run1()


if __name__ == '__main__':
    unittest.main()

