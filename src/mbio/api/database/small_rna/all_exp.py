# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from mbio.api.database.small_rna.api_base import ApiBase
import json
import pandas as pd
import datetime
from bson.objectid import ObjectId
import os
import numpy as np
from collections import OrderedDict
import math
from scipy import stats
import glob
import re
import unittest

class AllExp(ApiBase):
    def __init__(self, bind_object):
        super(AllExp, self).__init__(bind_object)

    ####################################################################################################

    @staticmethod
    def process_exp_matrix(exp_matrix, log_base=None, group_dict=None):

        # check exp_matrix type
        if type(exp_matrix) == str or type(exp_matrix) == bytes or isinstance(exp_matrix, unicode):
            all_exp_pd = pd.read_table(exp_matrix, index_col=0, header=0)
            all_exp_pd.index.name = 'seq_id'
        elif isinstance(exp_matrix, pd.DataFrame):
            self._bind_object.logger.debug('{} is assumed to be a pandas DataFrame Object'.format(exp_matrix))
            all_exp_pd = exp_matrix
            all_exp_pd.index.name = 'seq_id'
        else:
            self._bind_object.set_error('type of incoming exp_matrix is {}, abord'.format(type(exp_matrix)))

        # insert group mean columns
        if group_dict is not None:
            group_exp = list()
            for g in group_dict:
                g_exp = all_exp_pd.loc[:, group_dict[g]].mean(axis=1)
                g_exp.name = g
                group_exp.append(g_exp)
            all_exp_pd = pd.concat(group_exp, axis=1)

        # logarithmic transformation
        # important change: log(all_exp_pd + 1) --> log(all_exp_pd)
        if log_base:
            if log_base == math.e:
                all_exp_pd = np.log(all_exp_pd)
            elif log_base == 2:
                all_exp_pd = np.log2(all_exp_pd)
            elif log_base == 10:
                all_exp_pd = np.log10(all_exp_pd)
            else:
                self._bind_object.set_error('log base of {} is not supported'.format(log_base))

        return all_exp_pd

    ####################################################################################################

    def add_exp(self, exp_matrix, project_sn='small_rna', task_id='small_rna',
                exp_type='tpm', is_novel=False, params=None, libtype=None, main_id=None):

        # open exp_matrix as pd.DataFrame
        try:
            all_exp_pd = self.process_exp_matrix(exp_matrix)
        except Exception as e:
            self._bind_object.set_error('{} can not be processed by pd.DataFrame with exception: {}'.format(exp_matrix, e))

        # check params type
        if params == None:
            self._bind_object.logger.debug('type of incoming params is {}, pass'.format(type(params)))
        elif type(params) == dict:
            self._bind_object.logger.debug('type of incoming params is {}, dumps it to str'.format(type(params)))
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        elif type(params) == str:
            self._bind_object.logger.debug('type of incoming params is {}, pass'.format(type(params)))
        else:
            self._bind_object.set_error('type of incoming params is {}, abord'.format(type(params)))

        # if main_id is None, create a new document in sg_exp and return main_id
        if main_id == None:
            time_now = datetime.datetime.now()
            name = 'Exp_{}_{}'.format(exp_type, time_now.strftime('%Y%m%d_%H%M%S'))
            created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
            desc = 'exp_{}_main_table'.format(exp_type)
            # insert one document to sg_exp, return _id as main_id
            main_info = {
                'project_sn': project_sn,
                'task_id': task_id,
                'name': name,
                'exp_type': exp_type,
                'created_ts': created_ts,
                'desc': desc,
                'params': params,
                'libtype': libtype,
            }
            main_id = self.create_db_table('sg_exp', [main_info])
            self.update_db_record('sg_exp', main_id, main_id=main_id)
            self._bind_object.logger.info('succeed in creating document in sg_exp')

        # update exp_matrix path
        if is_novel:
            self.update_db_record('sg_exp', main_id, novel_matrix=exp_matrix)
        else:
            self.update_db_record('sg_exp', main_id, known_matrix=exp_matrix)

        # prepare record_dict_list and insert many documents to sg_exp_detail
        self.update_db_record('sg_exp', main_id, status='start')
        all_exp_pd['exp_id'] = ObjectId(main_id)
        all_exp_pd['is_novel'] = is_novel
        all_exp_pd.reset_index(level=0, inplace=True)
        row_dict_list = all_exp_pd.to_dict('records')
        # tranform all int to float
        for idx in range(len(row_dict_list)):
            for key in row_dict_list[idx].keys():
                if isinstance(row_dict_list[idx][key], int) and key != 'is_novel':
                    row_dict_list[idx][key] = float(row_dict_list[idx][key])
        self.create_db_table('sg_exp_detail', row_dict_list)
        self._bind_object.logger.info('succeed in creating documents in sg_exp_detail')
        self.update_db_record('sg_exp', main_id, status='end')

        # return main_id of sg_exp
        return main_id

    ####################################################################################################

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
            collection = self.db['sg_exp_pca_circ_detail']
            collection.insert_many(insert_data[1:])
        except Exception as e:
            self.bind_object.set_error("导入表格sg_exp_pca_circ_detail信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入sg_exp_pca_circ_detail表格成功")

    def add_distribution(self, exp_matrix, group_dict, project_sn='small_rna', task_id='small_rna',
                         exp_id=None, seq_type='all', exp_type='tpm', params=None, main_id=None):

        # open exp_matrix as pd.DataFrame, all_exp_pd no use here
        try:
            all_exp_pd = self.process_exp_matrix(exp_matrix)
        except Exception as e:
            self._bind_object.set_error('{} can not be processed by pd.DataFrame with exception: {}'.format(exp_matrix, e))

        self._bind_object.logger.debug('incoming main_id in add_distribution is {}'.format(main_id))

        # creat a new document if provide no main_id
        if main_id is None:

            # prepare main_info value
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            time_now = datetime.datetime.now()
            name = 'ExpDistribution_{}_{}_{}'.format(seq_type.lower(), exp_type, time_now.strftime('%Y%m%d_%H%M%S'))
            created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
            desc = 'stack and density and box and volin plot main table'

            # insert one document to sg_exp_graph, return _id as main_id
            main_info = {
                # essential keys
                'task_id': task_id,
                'project_sn': project_sn,
                'created_ts': created_ts,
                'name': name,
                'desc': desc,
                'params': params,
                # alterable keys
                'exp_id': ObjectId(exp_id),
                'seq_type': str(seq_type),
                # status of process
                'status': 'start'
            }
            main_id = self.create_db_table('sg_exp_graph', [main_info])
            self._bind_object.logger.info('succeed in creating documents in sg_exp_graph')
        else:
            self._bind_object.logger.info('skip creating documents in sg_exp_graph')

        # add graph info in detail
        main_id = ObjectId(main_id)
        self._bind_object.logger.info('start to add exp_graph_detail info to mongo')
        self._bind_object.logger.debug('incoming exp_matrix in add_distribution is {}'.format(exp_matrix))
        self._bind_object.logger.debug('type of incoming exp_matrix is {}'.format(type(exp_matrix)))
        self.add_stack(exp_matrix, group_dict, graph_id=main_id)
        self.add_box(exp_matrix, group_dict, graph_id=main_id)
        self.add_density(exp_matrix, group_dict, graph_id=main_id)
        self.add_volin(exp_matrix, group_dict, graph_id=main_id)
        self.update_db_record('sg_exp_graph', main_id, status='end', main_id=main_id)

        # return main_id of sg_exp_graph
        return main_id

    def add_stack(self, exp_matrix, group_dict, graph_id):
        '''
        dump exp data into sg_exp_graph_stack
        :param exp_matrix: expression matrix path or express matrix in pd.DataFrame format
        :param group_dict: group info dict
        :param graph_id: main table id
        '''

        # insert many documents to sg_exp_graph_volin with same graph_id and sample type
        all_exp_pd = self.process_exp_matrix(exp_matrix)
        row_dict_list = self.get_stack(all_exp_pd)
        self.create_db_table('sg_exp_graph_stack', row_dict_list, tag_dict={'graph_id': graph_id, 'type': 'sample'})

        # insert many documents to sg_exp_graph_volin with same graph_id and group type
        all_exp_pd = self.process_exp_matrix(exp_matrix, group_dict=group_dict)
        row_dict_list = self.get_stack(all_exp_pd)
        self.create_db_table('sg_exp_graph_stack', row_dict_list, tag_dict={'graph_id': graph_id, 'type': 'group'})
        self._bind_object.logger.info('succeed in creating documents in sg_exp_graph_volin')

    @staticmethod
    def get_stack(all_exp_pd):
        '''
        get stack plot info for each column of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        '''
        row_dict_list = list()
        for col in all_exp_pd.columns:
            milist = all_exp_pd.sort_values(by=col, ascending=False)[col]
            # stack_col = [{i: float(milist[i])/sum(milist)} for i in milist.index]
            # miRNA命名中包含特殊字符时不能作为mongo的键值
            stack_col = [{i.replace('.', '_'): milist[i]} for i in milist.index]
            row_dict_list.append({'name': col, 'stack_col': stack_col})
        return row_dict_list

    def add_box(self, exp_matrix, group_dict, graph_id):
        '''
        dump exp data into sg_exp_graph_box, only use log10 transformed values for drawing
        :param exp_matrix: expression matrix path or express matrix in pd.DataFrame format
        :param group_dict: group info dict
        :param graph_id: main table id
        :return: None
        '''

        # insert many documents to sg_exp_graph_box with same graph_id and sample type
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        stat_dict_list = self.get_box(all_exp_pd)
        self.create_db_table('sg_exp_graph_box', stat_dict_list, tag_dict={'graph_id': graph_id, 'type': 'sample'})

        # insert many documents to sg_exp_graph_box with same graph_id and group type
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        stat_dict_list = self.get_box(all_exp_pd)
        self.create_db_table('sg_exp_graph_box', stat_dict_list, tag_dict={'graph_id': graph_id, 'type': 'group'})
        self._bind_object.logger.info('succeed in creating documents in sg_exp_graph_box')

    @staticmethod
    def get_box(all_exp_pd):
        '''
        get box plot info for each column of the input pandas DataFrame
        :param all_exp_pd: pd.DataFrame
        :return: a list with dict as element
        '''
        stat_dict_list = list()
        target_columns = all_exp_pd.columns
        for each in target_columns:
            exp_pd = all_exp_pd[each]
            exp_pd = exp_pd[exp_pd != -np.inf]
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

    def add_density(self, exp_matrix, group_dict, graph_id):
        '''
        dump exp data into sg_exp_graph_density, only use log10 transformed values for drawing
        :param exp_matrix: expression matrix path or express matrix in pd.DataFrame format
        :param group_dict: group info dict
        :param graph_id: main table id
        '''

        # insert many documents to sg_exp_graph_density with same graph_id and sample type
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        records = self.get_density(all_exp_pd)
        self.create_db_table('sg_exp_graph_density', records, tag_dict={'graph_id': graph_id, 'type': 'sample'})

        # insert many documents to sg_exp_graph_density with same graph_id and group type
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        records = self.get_density(all_exp_pd)
        self.create_db_table('sg_exp_graph_density', records, tag_dict={'graph_id': graph_id, 'type': 'group'})
        self._bind_object.logger.info('succeed in creating documents in sg_exp_graph_density')

    @staticmethod
    def get_density(all_exp_pd):
        '''
        sampling 1000 density point for each columns of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        '''
        records = list()
        target_columns = all_exp_pd.columns
        for sample in target_columns:
            exp = all_exp_pd[sample]
            exp = exp[exp != -np.inf]
            density_func = stats.gaussian_kde(exp)
            # min_exp, max_exp = exp.min(), exp.max()
            # important change: min_exp = exp.min() --> min_exp = -2.0
            min_exp, max_exp = -2.0, exp.max()
            x_data = np.linspace(min_exp, max_exp, num=1000, endpoint=False)
            y_data = density_func(x_data)
            point_dict_list = pd.DataFrame({'lgexp': x_data, 'density': y_data}).to_dict('records')
            records.append({'sample': sample, 'data': point_dict_list})
        return records

    def add_volin(self, exp_matrix, group_dict, graph_id):
        '''
        dump exp data into sg_exp_graph_volin, only use log10 transformed values for drawing
        :param exp_matrix: expression matrix path or express matrix in pd.DataFrame format
        :param group_dict: group info dict
        :param graph_id: main table id
        '''

        # insert many documents to sg_exp_graph_volin with same graph_id and sample type
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10)
        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > -np.inf]
        all_exp_pd = all_exp_pd.reset_index(level=0)
        row_dict_list = all_exp_pd.to_dict('records')
        # row_dict_list = self.get_volin(all_exp_pd)
        self.create_db_table('sg_exp_graph_volin', row_dict_list, tag_dict={'graph_id': graph_id, 'type': 'sample'})

        # insert many documents to sg_exp_graph_volin with same graph_id and group type
        all_exp_pd = self.process_exp_matrix(exp_matrix, log_base=10, group_dict=group_dict)
        all_exp_pd[all_exp_pd == -np.inf] = 0
        all_exp_pd = all_exp_pd.reset_index(level=0)
        row_dict_list = all_exp_pd.to_dict('records')
        # row_dict_list = self.get_volin(all_exp_pd)
        self.create_db_table('sg_exp_graph_volin', row_dict_list, tag_dict={'graph_id': graph_id, 'type': 'group'})
        self._bind_object.logger.info('succeed in creating documents in sg_exp_graph_volin')

    @staticmethod
    def get_volin(all_exp_pd):
        '''
        sampling 70% point for each columns of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        '''
        all_exp_pd = all_exp_pd.sample(frac=0.7)
        all_exp_pd.reset_index(level=0, inplace=True)
        row_dict_list = all_exp_pd.to_dict('records')
        return row_dict_list

    ####################################################################################################

    def add_exp_venn(self, venn_graph, main_id=None, project_sn='small_rna', task_id='small_rna',
                     exp_level='A', method='miRDeep2', exp_type='tpm', exp_id=None, params=None, detail=True):
        '''
        create sg_exp_venn if main_id is not provide
        create sg_exp_venn_detail if assign detail as True by reading venn_graph
        :param venn_graph: venn_graph.xls resulted from exp_venn tool
        :param main_id: new document in sg_exp_venn would not been created when provide
        :param project_sn: project id
        :param task_id: task id
        :param exp_level: K or N or A
        :param method: exp quant method
        :param exp_type: tpm or count
        :param exp_id: main_id of related document in sg_exp
        :param params: string of parameters dict, designed for judgement of whether the task is repeated
        :return: main_id of sg_exp_venn
        '''

        # check whether maid_id is provided or not
        # if main_id is None, create a new document in sg_exp_venn and return main_id
        if main_id is None:
            time_now = datetime.datetime.now()
            name = 'ExpVenn_{}_{}_{}_{}'.format(exp_level, method, exp_type.upper(),
                                                time_now.strftime('%Y%m%d_%H%M%S'))
            created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
            if isinstance(params, dict):
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            desc = 'Expression venn main table built at api'
            main_info = {
                'project_sn': project_sn,
                'task_id': task_id,
                'name': name,
                'created_ts': created_ts,
                'exp_id': exp_id,
                'params': params,
                'desc': desc,
                'status': 'start'
            }
            main_id = self.create_db_table('sg_exp_venn', [main_info])

        # check whether detail is assigned or not
        # if detail is True, update sg_exp_venn after create a new document in sg_exp_venn_detail
        if detail:
            main_id = ObjectId(main_id)
            graph_pd = pd.read_table(venn_graph, header=0, sep='\t', keep_default_na = False)
            graph_pd.columns = ["name", "seqs"]
            detail_dict_list = graph_pd.to_dict('records')
            self.create_db_table('sg_exp_venn_detail', detail_dict_list, tag_dict={'venn_id': main_id})
            self.update_db_record('sg_exp_venn', record_id=main_id, status='end', main_id=main_id)
        return main_id

    ####################################################################################################

    def add_exp_corr(self, corr_output_dir, main_id=None, project_sn='small_rna', task_id='small_rna',
                     exp_level='A', method='miRDeep2', exp_type='tpm', exp_id=None, params=None, detail=True):
        '''
        create sg_exp_corr if main_id is not provide
        create sg_exp_corr_detail if assign detail as True by reading corr_output_dir
        :param corr_output_dir: output_dir resulted from exp_corr tool
        :param main_id: new document in sg_exp_corr would not been created when provide
        :param project_sn: project id
        :param task_id: task id
        :param exp_level: K or N or A
        :param method: exp quant method
        :param exp_type: tpm or count
        :param exp_id: main_id of related document in sg_exp
        :param params: string of parameters dict, designed for judgement of whether the task is repeated
        :return: main_id of sg_exp_corr
        '''

        # check whether maid_id is provided or not
        # if main_id is None, create a new document in sg_exp_corr and return main_id
        if main_id is None:
            time_now = datetime.datetime.now()
            name = 'ExpCorr_{}_{}_{}_{}'.format(exp_level, method, exp_type.upper(),
                                                time_now.strftime('%Y%m%d_%H%M%S'))
            created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
            if isinstance(params, dict):
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            desc = 'Expression corr main table built at api'
            main_info = {
                'task_id': task_id,
                'project_sn': project_sn,
                'created_ts': created_ts,
                'name': name,
                'desc': desc,
                'params': params,
                'exp_id': ObjectId(exp_id),
                'status': 'start'
            }
            main_id = self.create_db_table('sg_exp_corr', [main_info])

        # check whether detail is assigned or not
        # check corr_output_dir, read sample.cluster_tree.txt and sample_correlation.xls
        if detail:
            results = os.listdir(corr_output_dir)
            # cluster
            samples = list()
            sample_tree = ''
            if 'sample.cluster_tree.txt' in results:
                spct_file = os.path.join(corr_output_dir, 'sample.cluster_tree.txt')
                with open(spct_file) as f:
                    sample_tree = f.readline().strip()
                    samples = f.readline().strip().split(';')
            else:
                self._bind_object.set_error('fail to find sample.cluster_tree.txt in {}'.format(corr_output_dir))
            # corr
            ordered_dict_list = list()
            if 'sample_correlation.xls' in results:
                spcr_file = os.path.join(corr_output_dir, 'sample_correlation.xls')
                corr_result = pd.read_table(spcr_file, index_col=0, header=0)
                corr_result = corr_result.loc[samples, :]
                corr_result.reset_index(inplace=True)
                row_dict_list = corr_result.to_dict('records')
                ordered_dict_list = self.order_row_dict_list(row_dict_list, ['sample'] + samples)
            else:
                self._bind_object.set_error('fail to find sample_correlation.xls in {}'.format(corr_output_dir))
            # if detail is True, update sg_exp_corr after create a new document in sg_exp_corr_detail
            main_id = ObjectId(main_id)
            self.create_db_table('sg_exp_corr_detail', ordered_dict_list, tag_dict={'corr_id': main_id})
            self.update_db_record('sg_exp_corr', record_id=main_id,
                                  clust_tree=sample_tree, samples=samples, status="end", main_id=main_id)
            return main_id, sample_tree, samples
        else:
            return main_id

    ####################################################################################################

    def add_exp_pca(self, pca_output_dir, main_id=None, project_sn='small_rna', task_id='small_rna',
                    exp_level='A', method='miRDeep2', exp_type='tpm', exp_id=None, params=None, detail=True):
        '''
        create sg_exp_pca if main_id is not provide
        create sg_exp_pca_detail if assign detail as True by reading pca_output_dir
        :param pca_output_dir: output_dir resulted from exp_pca tool
        :param main_id: new document in sg_exp_pca would not been created when provide
        :param project_sn: project id
        :param task_id: task id
        :param exp_level: K or N or A
        :param method: exp quant method
        :param exp_type: tpm or count
        :param exp_id: main_id of related document in sg_exp
        :param params: string of parameters dict, designed for judgement of whether the task is repeated
        :return: main_id of sg_exp_pca
        '''

        # check whether maid_id is provided or not
        # if main_id is None, create a new document in sg_exp_pca and return main_id
        if main_id is None:
            time_now = datetime.datetime.now()
            name = 'ExpPca_{}_{}_{}_{}'.format(exp_level, method, exp_type.upper(),
                                                time_now.strftime('%Y%m%d_%H%M%S'))
            created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
            if isinstance(params, dict):
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            desc = 'Expression pca main table built at api'
            main_info = {
                'task_id': task_id,
                'project_sn': project_sn,
                'created_ts': created_ts,
                'name': name,
                'desc': desc,
                'params': params,
                'exp_id': exp_id,
                'status': 'start'
            }
            main_id = self.create_db_table('sg_exp_pca', [main_info])

        # check whether detail is assigned or not
        # check pca_output_dir, read PCA.xls and Explained_variance_ratio.xls
        if detail:
            main_id = ObjectId(main_id)
            pca_file = os.path.join(pca_output_dir, 'PCA.xls')
            result = pd.read_table(pca_file, header=0)
            row_dict_list = result.to_dict('records')
            self.create_db_table('sg_exp_pca_detail', row_dict_list, tag_dict={'pca_id': main_id})
            # if detail is True, update sg_exp_pca after create a new document in sg_exp_pca_detail
            evr_file = os.path.join(pca_output_dir, 'Explained_variance_ratio.xls')
            t = pd.read_table(evr_file, header=None)
            pc_ratio_dict = OrderedDict(zip(t[0], t[1]))
            self.update_db_record('sg_exp_pca', record_id=main_id,
                                  ratio_dict=pc_ratio_dict, status="end", main_id=main_id)
        return main_id

    ####################################################################################################

    def add_diffexp(self, diff_output, main_id=None, exp_id=None, group_dict=None, group_id=None,
                    project_sn='small_rna', task_id='small_rna', params=None,
                    exp_level='M', method='miRDeep2', exp_type='tpm',diff_method='DESeq2',
                    pvalue_padjust='padjust', create_geneset=True):
        '''
        add differential analysis result to database
        :param diff_output: diffexp result dir
        :param exp_id: exp table id from POST, will be used for getting expression value
        :param group_dict: group info dict. {group:[s1,s2,], ...}
        :param exp_level: expression level
        :param pvalue_padjust: pvalue or padjust, for significant judgement.
        :param method: method for expression quant
        :param method: tpm
        :param diff_method: differential analysis method
        :param project_sn: project id
        :param task_id: task id
        :param main_id: new document in sg_diff would not been created when provide
        :param group_id: main_id of related document in sg_specimen_group
        :param create_geneset: whether create DE geneset or not
        :param params: from POST, designed for judgement of whether the task is repeated.
        :return: main_id of sg_diff
        '''

        summary = glob.glob(os.path.join(diff_output, '*_diff_summary.xls'))[0]
        diff_files = glob.glob(os.path.join(diff_output, '*_vs_*.*.xls'))
        if not diff_files:
            self._bind_object.set_error('no target file found in {}'.format(diff_output))

        # check whether maid_id is provided or not
        # if main_id is None, create a new document in sg_diff and return main_id
        if main_id is None:
            time_now = datetime.datetime.now()
            name = 'Diffexp_{}_{}'.format(diff_method, time_now.strftime('%Y%m%d_%H%M%S'))
            created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
            if isinstance(params, dict):
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            desc = 'Differential expression main table built at api'
            main_info = {
                # essential keys
                'task_id': task_id,
                'project_sn': project_sn,
                'created_ts': created_ts,
                'name': name,
                'desc': desc,
                'params': params,
                # alterable keys
                'exp_id': exp_id,
                # status of process
                'status': 'start'
            }
            main_id = self.create_db_table('sg_diff', [main_info])

        main_id = ObjectId(main_id)
        cmp_list = list()
        cmp_detail_dict = dict()
        diff_dict_list = list()
        volcano_dict_list = list()
        sig_status = dict()
        scatter_dict_list = list()

        # add key 'is_novel' to documents in sg_diff_detail
        seq_type_pd = self.get_known_seq_list(exp_id=exp_id)

        for each in diff_files:
            if each.endswith('.annot.xls'):
                continue
            if each.endswith('.normalize.xls'):
                continue
            if each.endswith('.sizeFactor.xls'):
                continue
            if each.endswith('.normFactor.xls'):
                continue
            diff_pd = pd.read_table(each, header=0, sep='\t')
            columns = diff_pd.columns
            log2fc_ind = list(columns).index('log2fc')
            need_cols = ['seq_id', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
            need_cols += [columns[log2fc_ind-2], columns[log2fc_ind-1]]
            samples = list()
            for x in columns:
                m = re.match(r'(.*)_count$', x)
                if m:
                    samples.append(m.groups()[0])
            fname = os.path.basename(each)
            ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(diff_method.lower()), fname).groups()
            cmp_combine = ctrl + '|' + test
            cmp_list.append(cmp_combine)
            cmp_detail_dict[cmp_combine] = samples
            cmp_pd = pd.DataFrame([cmp_combine]*diff_pd.shape[0], columns=['compare'])
            tmp_pd = pd.concat([diff_pd.loc[:, need_cols], cmp_pd], axis=1)
            tmp_pd.columns = list(tmp_pd.columns[:-3]) + ['group1', 'group2', 'compare']

            # add key 'is_novel' to documents in sg_diff_detail
            tmp_pd = pd.merge(tmp_pd, seq_type_pd)

            diff_dict_list = tmp_pd.to_dict('records')
            # add detail info
            self.create_db_table('sg_diff_detail', diff_dict_list, tag_dict={'diff_id': main_id})

            # get volcano data
            status_list, stat_cutoff = self.get_volcano_status_cutoff(diff_pd, pvalue_padjust)
            sig_status[cmp_combine] = status_list
            volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', pvalue_padjust, 'significant', 'regulate']]
            # important change: remove 1og10pvalue adjustment
            # bool_ind = volcano_pd[pvalue_padjust] <= 0
            # volcano_pd.loc[bool_ind, pvalue_padjust] = 1e-200
            volcano_pd[pvalue_padjust] = -np.log10(volcano_pd[pvalue_padjust])
            lg_p_series = volcano_pd[pvalue_padjust]
            # transfer inf in log10pvalue to max rational number
            lg_p_max = max(lg_p_series[lg_p_series != np.inf])
            lg_p_series[lg_p_series == np.inf] = lg_p_max
            volcano_pd[pvalue_padjust] = lg_p_series
            volcano_pd.dropna(inplace=True)
            volcano_pd.columns = ['seq_id', 'log2fc', 'log10pvalue', 'significant', 'regulate']
            bool_idx = volcano_pd['log10pvalue'] <= 0
            volcano_pd.loc[bool_idx, 'log10pvalue'] = 0.0
            # important change: remove 1og10pvalue adjustment
            # bool_ind = volcano_pd['log10pvalue'] > stat_cutoff
            # volcano_pd.loc[bool_ind, 'log10pvalue'] = stat_cutoff
            volcano_pd = pd.concat([volcano_pd, cmp_pd], axis=1)
            volcano_pd_nosig = volcano_pd[volcano_pd['significant'] == 'no']
            # important change: remove random selection
            # random select 10000 not sig diff genes for plotting
            # if volcano_pd_nosig.shape[0] > 10000:
            #     volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
            #     volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.5)
            #     if volcano_pd_nosig.shape[0] > 15000:
            #         volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
            #     volcano_pd = pd.concat([volcano_pd_sig, volcano_pd_nosig], axis=0)
            volcano_dict_list += volcano_pd.to_dict('records')

            # get scatter data
            scatter_pd = tmp_pd.loc[:, ['seq_id', 'group1', 'group2', 'compare', 'significant', 'regulate']]
            scatter_pd.set_index('seq_id', inplace=True)
            scatter_pd = scatter_pd.loc[volcano_pd['seq_id'], :].reset_index()
            # important change: log(group + 1) --> log(group) and set min cutoff as -2.0
            # scatter_pd['group1'] = (scatter_pd['group1']+1).apply(np.log10)
            # scatter_pd['group2'] = (scatter_pd['group2']+1).apply(np.log10)
            for group in ['group1', 'group2']:
                tmp_col = scatter_pd[group].apply(np.log10)
                tmp_col[tmp_col <= -2] = -2.0
                scatter_pd[group] = tmp_col
            scatter_dict_list += scatter_pd.to_dict('records')

            # get significant diff gene set and add to database
            if create_geneset:
                sig_seqs = list(diff_pd['seq_id'][diff_pd['significant'] == 'yes'])
                geneset_main_info = dict(
                    project_sn=project_sn,
                    task_id=task_id,
                    name='{}_vs_{}'.format(ctrl, test),
                    type=exp_level,
                    desc='Differential expressed gene set',
                    group_id=group_id,
                    gene_length=len(sig_seqs),
                    is_use=1
                )
                genet_detail_info = [{"seq_list": sig_seqs, }]
                if len(sig_seqs) != 0:
                    self.add_set(geneset_main_info, genet_detail_info)
        else:
            self._bind_object.logger.info('loop end of diff_files')

        # add volcano detail
        self.create_db_table('sg_diff_volcano', volcano_dict_list, tag_dict={'diff_id': main_id})

        # add scatter detail
        self.create_db_table('sg_diff_scatter', scatter_dict_list, tag_dict={'diff_id': main_id})

        # assign first and second rows as header (two level)
        summary_pd = pd.read_table(summary, header=[0, 1])
        levels = summary_pd.columns.levels
        labels = summary_pd.columns.labels

        # add summary detail
        summary_pd.columns = levels[0][labels[0]]

        # add key 'is_novel' to documents in sg_diff_summary
        summary_pd = pd.merge(summary_pd, seq_type_pd)

        summary_dict_list = summary_pd.to_dict('records')
        self.create_db_table('sg_diff_summary', summary_dict_list, tag_dict={'diff_id': main_id})

        # update status
        self.update_db_record('sg_diff', record_id=main_id, cmp_list=cmp_list, cmp_detail=cmp_detail_dict,
                              sig_status=sig_status, status="end", main_id=main_id)
        return main_id

    @staticmethod
    def get_volcano_status_cutoff(diff_table, pvalue_padjust):
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
        log10_pvalue_list = list(log10_sig_pvalues[log10_sig_pvalues > 0])

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
            if len(log10_pvalue_list) > 1:
                log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.8)]
            else:
                log10_pvalue_cutoff = 200
        return sig_status, log10_pvalue_cutoff

    def get_known_seq_list(self, exp_id):
        conn = self.db['sg_exp_detail']
        results = conn.find({'exp_id': ObjectId(exp_id)}, {'seq_id': 1, 'is_novel': 1})
        data = {'seq_id': list(), 'is_novel': list()}
        for result in results:
            data['seq_id'].append(result['seq_id'])
            data['is_novel'].append(result['is_novel'])
        ret = pd.DataFrame(data)
        return ret

    def add_set(self, main_info, detail_info):
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        main_info.update(dict(status='start', created_ts=created_ts))
        main_id = self.create_db_table('sg_geneset', [main_info])
        self.create_db_table('sg_geneset_detail', detail_info, tag_dict={'geneset_id': main_id})
        task_id = main_info['task_id']
        self.update_db_record('sg_geneset', record_id=main_id, status='end', is_use=0, main_id=main_id)
        return main_id

    ####################################################################################################

    def run(self):
        project_sn = "55617_5f3e322484a2a"
        task_id = "majorbio_279548"
        known_count_xls = "/mnt/ilustre/users/isanger/workspace/20200825/Smallrna_majorbio_279548/Srna/output/known_mirna_count.xls"
        novel_count_xls = "/mnt/ilustre/users/isanger/workspace/20200825/Smallrna_majorbio_279548/Srna/output/novel_mirna_count.xls"
        known_tpm_xls = "/mnt/ilustre/users/isanger/workspace/20200825/Smallrna_majorbio_279548/Srna/output/known_mirna_norm.xls"
        novel_tpm_xls = "/mnt/ilustre/users/isanger/workspace/20200825/Smallrna_majorbio_279548/Srna/output/novel_mirna_norm.xls"
        count_main_id = self.add_exp(
            exp_matrix=known_count_xls,
            project_sn=project_sn,
            task_id=task_id,
            exp_type='count',
            is_novel=False
        )
        count_exp_id = self.add_exp(
            exp_matrix=novel_count_xls,
            project_sn=project_sn,
            task_id=task_id,
            exp_type='count',
            is_novel=True,
            main_id=str(count_main_id)
        )

class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test(test):
        from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "majorbio_279548",
            "project_sn": "55617_5f3e322484a2a",
            "type": "workflow",
            "name": "majorbio_279548",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = SmallRnaTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.api = wf.api.api("small_rna.all_exp")
        wf.api.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)