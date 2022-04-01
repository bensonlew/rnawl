# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId

import json
import os
import unittest
from collections import OrderedDict
from biocluster.config import Config
import pandas as pd
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from mbio.api.database.whole_transcriptome.api_base import ApiBase



class Batch(ApiBase):
    def __init__(self, bind_object):
        super(Batch, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'


    def add_exp_batch(self, count_batch, exp_all, exp_level, main_exp_id, exp_id, sg_exp_batch_id, batch_version, task_id,params_json, category, library):
        if type(main_exp_id) == str or type(main_exp_id) == bytes or type(main_exp_id) == unicode:
            main_exp_id = ObjectId(main_exp_id)
        if type(exp_id) == str or type(exp_id) == bytes or type(exp_id) == unicode:
            exp_id = ObjectId(exp_id)
        connect = self.db['exp_detail']
        now_batch_version = connect.find_one({'exp_id': exp_id})['batch_version']
        df = pd.read_table(count_batch, header=0, index_col=False, sep='\t')
        df2 = pd.read_table(exp_all, sep='\t', header=0)
        all_columns = df2.columns.tolist()
        if exp_level == 'T':
            exp_all_matrix = df2.set_index('transcript_id')
            df.rename(columns={'seq_id': 'transcript_id'}, inplace=True)
            exp_matrix = df.set_index('transcript_id')
        if exp_level == 'G':
            exp_all_matrix = df2.set_index('gene_id')
            df.rename(columns={'seq_id': 'gene_id'}, inplace=True)
            exp_matrix = df.set_index('gene_id')

        columns_dict = dict()
        for i in exp_matrix.columns.tolist():
            columns_dict[i] = '{}_batch'.format(i)
        if batch_version == 0:
            for i in exp_matrix.columns.tolist():
                exp_all_matrix['{}_batch'.format(i)] = exp_all_matrix[i]
        if batch_version > 0:
            pass
        exp_matrix.rename(columns=columns_dict, inplace=True)
        exp_all_matrix.update(exp_matrix)
        exp_final = exp_all_matrix.reset_index()



        # if batch_version == 0:
        #     db_matrix = pd.concat([exp_all_matrix, exp_matrix], axis=1)
        #     if exp_level == 'G':
        #         db_matrix.index.name = 'gene_id'
        #     else:
        #         db_matrix.index.name = 'transcript_id'
        #     exp_final = db_matrix.reset_index()
        # else:
        #     exp_all_matrix.update(exp_matrix)
        #     exp_final = exp_all_matrix.reset_index()
        # exp_all_matrix = exp_all_matrix.reset_index()
        # exp_matrix = exp_matrix.reset_index()
        # if exp_level == 'T':
        #     db_matrix = pd.merge(exp_matrix,exp_all_matrix, how='outer', on='transcript_id')
        # if exp_level == 'G':
        #     db_matrix = pd.merge(exp_matrix, exp_all_matrix, how='outer', on='gene_id')
        # pd.merge()



        # _id = exp_final['_id']
        _exp_id = exp_final['exp_id']
        _exp_id_new = [ObjectId(x) for x in _exp_id]
        exp_final['exp_id'] = _exp_id_new
        # _id_new = [ObjectId(x) for x in _id]
        batch_version_new = batch_version + 1
        exp_final['batch_version'] = batch_version_new
        if batch_version_new == now_batch_version:
            raise OptionError("不能插入batch_version相同的表达量详情表", code = "35002901")
        exp_final = exp_final.fillna("")
        final_columns = exp_final.columns.tolist()
        # print final_columns, all_columns
        # if all_columns == final_columns:
        #     pass
        # else:
        #     raise OptionError("插入的列数与原始的列数不同", code = "35002901")
        detail = exp_final.to_dict('r')
        # connect = self.db['sg_exp_detail']
        # for i in detail:
        #     transcript_id = i['seq_id']
        #     connect.update({'exp_id': exp_id, 'transcript_id': transcript_id}, {'$set': i})

        self.create_db_table('exp_detail', detail)
        self.remove_db_record('exp_detail', exp_id=exp_id, batch_version=batch_version)

        # db = Config().get_mongo_client(mtype=self._project_type)[Config().get_mongo_dbname(self._project_type)]
        # connect = db['sg_exp']
        # record = connect.find_one({'main_id': main_exp_id, 'is_rmbe': 'true', 'status': 'end'})
        # record.remove()
        db = Config().get_mongo_client(mtype=self._project_type)[Config().get_mongo_dbname(self._project_type)]
        try:
            record = db['exp'].find_one({'is_rmbe': 'true', 'batch_main_id': exp_id, "status": 'end'})
            category_old = record['category']
            category_new = category_old + ';' + category
            library_old = record['library']
            library_new = library_old + ';' + library
        except:
            category_new = category
            library_new = library



        try:
            self.remove_db_record('exp', is_rmbe='true', batch_main_id=exp_id, status='end')
        except:
            pass
        self.update_db_record('exp', exp_id, batch_version=batch_version_new)
        self.update_db_record('exp', main_exp_id, status='end', batch_version=batch_version_new, category=category_new, library=library_new)
        # self.update_db_record('sg_exp', main_exp_id, main_id=exp_id)
        self.update_db_record('exp', main_exp_id, batch_main_id=exp_id)
        # self.update_db_record('sg_exp_batch', sg_exp_batch_id, insert_dict={'main_id': sg_exp_batch_id})

        # self.create_db_table('sg_exp_detail', exp_matrix.to_dict('r'))
        # self.update_db_record('sg_exp_batch', sg_exp_batch_id, insert_dict={'main_id': sg_exp_batch_id})


    def add_exp_batch_corr(self, corr_output_dir, batch_version, task_id, library, level, main_id=None, record_id=None):
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
        batch_version_new = batch_version + 1
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
        # if not isinstance(record_id, ObjectId):
        #     record_id = ObjectId(record_id)
        self.update_db_record('exp_batch', main_id, clust_tree=sample_tree, samples=samples, )
        tag_dict = dict(batch_id=main_id, batch_version=batch_version_new)
        self.create_db_table('exp_batch_corr_detail', ordered_dict_list, tag_dict=tag_dict)
        self.update_db_record('exp_batch', main_id, batch_version=batch_version_new)

        db = Config().get_mongo_client(mtype=self._project_type)[Config().get_mongo_dbname(self._project_type)]
        connect = db['exp_batch']
        # record = connect.find_one({'task_id': task_id, 'batch_version': batch_version, 'library': library})
        record = connect.find_one({'task_id': task_id, 'level': level, 'library': library, 'status': 'end'})
        if record:
            corr_id = record['main_id']
            self.remove_db_record('exp_batch', main_id=corr_id)
            self.remove_db_record('exp_batch_corr_detail', batch_id=corr_id)
        else:
            pass
        self.update_db_record('exp_batch', main_id, status="end")

        return main_id, sample_tree, samples

    def add_exp_batch_pca(self, exp_output_dir, batch_version, task_id, main_id, record_id, library, level):
        batch_version_new = batch_version+1
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
        self.update_db_record('exp_batch', main_id, ratio_dict=pc_ratio_dict)
        # insert detail
        tag_dict = dict(batch_id=main_id, batch_version=batch_version_new)
        self.create_db_table('exp_batch_pca_detail', row_dict_list, tag_dict=tag_dict)
        self.update_db_record('exp_batch', main_id, batch_version=batch_version_new)
        db = Config().get_mongo_client(mtype=self._project_type)[Config().get_mongo_dbname(self._project_type)]
        connect = db['exp_batch']
        # record = connect.find_one({'task_id': task_id, 'batch_version': batch_version, 'library': library})
        record = connect.find_one({'task_id': task_id, 'level': level, 'library': library, 'status': 'end'})
        if record:
            pca_id = record['main_id']
            # self.remove_db_record('sg_exp_batch',main_id=pca_id)
            self.remove_db_record('exp_batch_pca_detail', batch_id=pca_id)
        else:
            pass
        return main_id

    def insert_ellipse_table(self, infile, batch_version, task_id, main_id, record_id, library, level):
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
                    tmp['batch_id'] = main_id
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
            collection = self.db['exp_batch_pca_circ_detail']
            collection.insert_many(insert_data[1:])
            db = Config().get_mongo_client(mtype=self._project_type)[Config().get_mongo_dbname(self._project_type)]
            connect = db['exp_batch']
            # record = connect.find_one({'task_id': task_id, 'batch_version': batch_version, 'library': library})
            record = connect.find_one({'task_id': task_id, 'level': level, 'library': library, 'status': 'end'})
            if record:
                exp_pca_ellipse_id = record['main_id']
                self.remove_db_record('exp_batch_corr_detail', batch_id=exp_pca_ellipse_id)
            else:
                pass
        except Exception as e:
            self.bind_object.set_error("导入表格sg_exp_batch_pca_circ_detail信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入sg_exp_batch_pca_circ_detail表格成功")