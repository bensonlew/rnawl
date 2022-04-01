# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.dia.api_base import ApiBase
import re
import os
import glob
import pandas as pd
from bson.objectid import ObjectId
from collections import OrderedDict
import numpy as np
import types


class Preprocess(ApiBase):
    def __init__(self, bind_object):
        super(Preprocess, self).__init__(bind_object)
        self._project_type = 'dia'

    #@report_check
    def add_preprocess_exp(self, exp_path, cv_path, cv_summary, searchdb, exp_type='origin', raw_exp_id=None,
                           nrmse_path=None, na_path=None, main_id=None, project_sn=None, task_id=None, params=None):
        if main_id is None:
            # prepare main table info
            if exp_type == 'raw':
                name = 'Raw'
            else:
                name = "ProteinTable_Origin"
            time_now = datetime.datetime.now()
            # name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='exp main table',
                params=params,
                status="start",
                version='v3',
                type=exp_type,
            )
            try:
                main_id = self.create_db_table('sg_express', [main_info])
                self.bind_object.logger.info("导入sg_express主表{}成功".format(main_id))
            except Exception as e:
                self.bind_object.logger.info("导入sg_express主表信息出错:%s" % (e,))
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        if exp_type == 'raw':
            self.add_preprocess_detail(main_id, exp_path, cv_path, cv_summary, searchdb=searchdb, raw_cv=True)
        else:
            self.add_preprocess_detail(main_id, exp_path, cv_path, cv_summary, searchdb=searchdb, nrmse_path=nrmse_path,
                                       raw_exp_id=raw_exp_id, na_path=na_path)
        return main_id

    def add_preprocess_detail(self, express_id, exp_path, cv_path, cv_summary, searchdb, nrmse_path=None,
                              raw_exp_id=None, na_path=None, interactive=None, raw_cv=False):
        if not isinstance(express_id, ObjectId):
            if isinstance(express_id, types.StringTypes):
                express_id = ObjectId(express_id)
            else:
                self.bind_object.set_error('disgenet_enrich_id必须为ObjectId对象或其对应的字符串!', )

        if os.path.exists(searchdb):
            searchdb_pd = pd.read_table(searchdb, index_col=False, header=0)
            searchdb_selected = pd.DataFrame(index=searchdb_pd.index)
            for col in ['Accession', 'Description']:
                try:
                    searchdb_selected[col] = searchdb_pd[col]
                except:
                    try:
                        searchdb_selected[col] = searchdb_pd['# ' + col]  # compatible to self.new in main workflow
                    except:
                        self.bind_object.logger.info('protein文件里缺少%s这一列，请检查' % col)
                        searchdb_selected[col] = '_'

        if os.path.exists(exp_path) and os.path.exists(cv_path) and os.path.exists(searchdb):
            exp_pd = pd.read_table(exp_path, index_col=False, header=0)
            samples = list(exp_pd.columns[1:])
            cv_pd = pd.read_table(cv_path, index_col=False, header=0)
            groups = list(cv_pd.columns[1:-1])
            final_pd = pd.merge(cv_pd, exp_pd, on="Accession")
            final_pd = pd.merge(final_pd, searchdb_selected, on="Accession")
            final_pd.rename(columns={final_pd.columns[0]: "accession_id"}, inplace=True)
            final_pd.rename(columns={"Description": 'description'}, inplace=True)
            final_pd["express_id"] = express_id
            exp_dict_list = final_pd.to_dict('records')
            try:
                self.create_db_table('sg_express_detail', exp_dict_list)
            except Exception as e:
                self.bind_object.logger.info("导入sg_express_detail:%s信息出错:%s" % (express_id,e,))
            else:
                self.bind_object.logger.info("导入sg_express_detail:%s成功" % (express_id,))
            self.add_cv_box(cv_matrix=cv_path, express_id=express_id, raw_cv=raw_cv)

        if os.path.exists(cv_summary):
            cv_sum_pd = pd.read_table(cv_summary, index_col=0, header=0)
            col = list(cv_sum_pd.columns)
            cv_data = list()
            for g in groups:
                f_dict = dict()
                p_dict = dict()
                c_dict = dict()
                cat_list = list()
                for i in col[:-1]:
                    term = i.split("_")[1]
                    category = i.split("_")[0].split("%")[0]
                    if term == "frequency":
                        frequency = cv_sum_pd.loc[[g], [i]]
                        f_dict[category] = frequency.iloc[0][0]
                        cat_list.append(category)
                    if term == "percent":
                        percent = cv_sum_pd.loc[[g], [i]]
                        p_dict[category] = percent.iloc[0][0]
                    if term == "cumulative":
                        cumulative = cv_sum_pd.loc[[g], [i]]
                        c_dict[category] = cumulative.iloc[0][0]
                mean = cv_sum_pd.loc[[g], [cv_sum_pd.columns[-1]]]
                mongo_data = {
                    "express_id": express_id,
                    "frequency": f_dict,
                    "percent_frequency": p_dict,
                    "cumulative": c_dict,
                    "category": cat_list,
                    "group": g.split('cv_')[1],
                    "mean_cv": mean.iloc[0][0]
                }
                cv_data.append(mongo_data)
            try:
                collection = self.db["sg_express_cv_graph"]
                collection.insert_many(cv_data)
            except Exception as e:
                self.bind_object.logger.info("导入sg_express_cv_graph:%s信息出错:%s" % (express_id, e,))
            else:
                self.bind_object.logger.info("导入sg_express_cv_graph:%s成功" % (express_id,))

        if interactive != 'yes' and na_path:
            na_pd = pd.read_table(na_path, index_col=False, header=0)
            na_pd = pd.merge(na_pd, searchdb_selected, on="Accession")
            na_pd.rename(columns={na_pd.columns[0]: "accession_id"}, inplace=True)
            na_pd.rename(columns={"Description": 'description'}, inplace=True)
            if type(raw_exp_id) == str or type(raw_exp_id) == bytes or type(raw_exp_id) == unicode:
                raw_exp_id = ObjectId(raw_exp_id)
            na_pd["express_id"] = raw_exp_id
            na_pd.reset_index(level=0, inplace=True)
            na_dict_list = na_pd.to_dict('records')
            try:
                self.create_db_table('sg_express_na', na_dict_list)
            except Exception as e:
                self.bind_object.logger.info("导入sg_express_na:%s信息出错:%s" % (express_id, e,))
            else:
                self.bind_object.logger.info("导入sg_express_na:%s成功" % (express_id,))

        if nrmse_path:
            fill_na = os.path.basename(nrmse_path).split("_")[0]
            nrmse_pd = pd.read_table(nrmse_path, index_col=False, header=0)
            nrmse_pd.rename(columns={nrmse_pd.columns[0]: "fill_na"}, inplace=True)
            fill_nrmse = nrmse_pd['fill_na'][0]
            zero_nrmse = nrmse_pd['zero'][0]
            self.update_db_record('sg_express', express_id, method=fill_na, fill_na=fill_nrmse, zero=zero_nrmse)
            self.bind_object.logger.info("更新nrmse信息成功，sg_express:%s" % (express_id,))

        group = [each.split('cv_')[1] for each in groups]
        self.update_db_record('sg_express', express_id, status="end", main_id=express_id,
                              samples=samples, groups=group)

    @staticmethod
    def get_cv_box(group_cv_pd):
        """
        get box plot info for each column of the input pandas DataFrame
        :param group_cv_pd: pandas DataFrame
        :return: a list with dict as element
        """
        stat_dict_list = list()
        target_columns = group_cv_pd.columns
        for each in target_columns:
            cv_pd = group_cv_pd[each]
            cv_pd = cv_pd[cv_pd != 0]
            summary = cv_pd.describe()
            summary.index = [u'count', u'mean', u'std', u'min', u'q1', u'median', u'q3', u'max']
            tmp_dict = summary.to_dict()
            lt25 = cv_pd[cv_pd <= tmp_dict['q1']].shape[0]
            lt50 = cv_pd[cv_pd <= tmp_dict['median']].shape[0]
            lt75 = cv_pd[cv_pd <= tmp_dict['q3']].shape[0]
            upper_whisker = tmp_dict['q3'] + 1.5*(tmp_dict['q3'] - tmp_dict['q1'])
            lower_whisker = tmp_dict['q1'] - 1.5*(tmp_dict['q3'] - tmp_dict['q1'])
            if lower_whisker < 0:
                lower_whisker = 0
            upper_outliers = list(cv_pd[cv_pd > upper_whisker])
            lower_outliers = list(cv_pd[cv_pd < lower_whisker])
            tmp_dict.update({
                'group': each.split('cv_')[1],
                'min-q1': lt25,
                'q1-median': lt50-lt25,
                'median-q3': lt75-lt50,
                'q3-max': cv_pd.shape[0]-lt75,
                'upper_whisker': upper_whisker,
                'lower_whisker': lower_whisker,
                'upper_outliers': upper_outliers,
                'lower_outliers': lower_outliers,
            })
            stat_dict_list.append(tmp_dict)
        return stat_dict_list

    def add_cv_box(self, cv_matrix, express_id, raw_cv):
        """
        dump exp data into database
        :param cv_matrix: cv stats matrix path or an cv matrix in pandas DataFrame format.
        :param express_id: main table id
        """
        # for group info
        # if type(cv_matrix) == str or type(cv_matrix) == bytes or isinstance(cv_matrix, unicode):
        #     group_cv_pd = pd.read_table(cv_matrix, index_col=0, header=0)
        # else:
        #     print(cv_matrix, 'is assumed to be a pandas DataFrame Object')
        #     group_cv_pd = cv_matrix
        # group_cv_pd.index.name = 'accession_id'
        # group_cv_pd = group_cv_pd[group_cv_pd.sum(axis=1) > 0.001]
        # stat_list = self.get_cv_box(group_cv_pd)
        if type(express_id) == str or type(express_id) == bytes or type(express_id) == unicode:
            express_id = ObjectId(express_id)
        if raw_cv:
            box_file = os.path.join(os.path.dirname(cv_matrix), 'boxdata_raw.xls')
        else:
            box_file = os.path.join(os.path.dirname(cv_matrix), 'boxdata.xls')
        if os.path.exists(box_file):
            box_pd = pd.read_table(box_file, index_col=False, header=0, keep_default_na=False)
            box_pd['upper_outliers'] = box_pd['upper_outliers'].apply(lambda x: map(eval, str(x).split(';')) if x else [])
            box_pd['lower_outliers'] = box_pd['lower_outliers'].apply(lambda x: map(eval, str(x).split(';')) if x else [])
            stat_list = box_pd.to_dict('records')
            self.create_db_table('sg_express_cv_box', stat_list, tag_dict=dict(express_id=express_id))
