# !/usr/bin/python
# -*- coding: utf-8 -*-
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


class Diffexp(ApiBase):
    def __init__(self, bind_object):
        super(Diffexp, self).__init__(bind_object)
        self._project_type = 'tool_lab'
    @staticmethod
    def _get_volcano_status_cutoff(diff_table, pvalue_padjust):
        sig_status = list()
        sig_mark = diff_table['significant']
        reg_list = diff_table['regulate']
        down_num = 0
        up_num = 0
        if 'no' in list(sig_mark):
            no_sig_num = sig_mark[sig_mark == "no"].shape[0]
            sig_status.append('nosig_' + str(no_sig_num))
        if 'yes' in list(sig_mark):
            reg_mark = reg_list[sig_mark == 'yes']
            if 'down' in list(reg_mark):
                down_num = reg_mark[reg_mark == 'down'].shape[0]
            sig_status.append('down_' + str(down_num))
            if 'up' in list(reg_mark):
                up_num = reg_mark[reg_mark == 'up'].shape[0]
            sig_status.append('up_' + str(up_num))

        sig_pvalues = diff_table[pvalue_padjust][diff_table['significant'] == "yes"]
        log10_sig_pvalues = -np.log10(sig_pvalues)
        log10_pvalue_list = sorted(list(log10_sig_pvalues[log10_sig_pvalues > 0]))

        if len(sig_pvalues) > 2000:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.85)]
        elif len(sig_pvalues) > 1000:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.90)]
        elif len(sig_pvalues) > 500:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.95)]
        elif len(sig_pvalues) > 250:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.99)]
        elif len(sig_pvalues) == 0:
            tmp = -np.log10(diff_table[pvalue_padjust])
            tmp_list = sorted(tmp[tmp > 0])
            if len(tmp_list) == 0:
                log10_pvalue_cutoff = 200
            else:
                log10_pvalue_cutoff = tmp_list[int(len(tmp_list) * 0.9)]
        else:
            # print(pvalue_padjust, diff_table, log10_pvalue_list)
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.8)]
        return sig_status, log10_pvalue_cutoff, down_num, up_num

    def _get_volcano_status_cutoff_noiseq(self, diff_table):
        sig_status = list()
        sig_mark = diff_table['significant']
        reg_list = diff_table['regulate']
        down_num = 0
        up_num = 0
        if 'no' in list(sig_mark):
            no_sig_num = sig_mark[sig_mark == "no"].shape[0]
            sig_status.append('nosig_' + str(no_sig_num))
        if 'yes' in list(sig_mark):
            reg_mark = reg_list[sig_mark == 'yes']
            if 'down' in list(reg_mark):
                down_num = reg_mark[reg_mark == 'down'].shape[0]
            sig_status.append('down_' + str(down_num))
            if 'up' in list(reg_mark):
                up_num = reg_mark[reg_mark == 'up'].shape[0]
            sig_status.append('up_' + str(up_num))

        return sig_status, down_num, up_num

    def add_diffexp(self, diff_output
                    , quant_method='RSEM', diff_method='DESeq2', main_id=None,
                    project_sn='denovo_rna_v2', task_id='denovo_rna_v2', params=None,
                    pvalue_padjust='padjust', create_geneset=True, nosig_info=None, s3_output=None, pdf_name=None):
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
            # self.create_db_table('{}_diff_detail'.format(diff_method), nosig_list, tag_dict={'diff_id': main_id, "seq_type": "no_sig"})
        summary = glob.glob(os.path.join(diff_output, '*diff_summary*.xls'))[0]
        # diff_method = os.path.basename(summary).split('_diff_summary.xls')[0]
        diff_files = glob.glob(os.path.join(diff_output, '*_vs_*.*.xls'))
        if not diff_files:
            self.bind_object.set_error('No target file found in %s', variables=(diff_output), code="53700204")
        cmp_list = list()
        # cmp_detail_dict = dict()
        # diff_dict_list = list()
        volcano_dict_list = list()
        # sig_status = dict()
        # scatter_dict_list = list()
        value_list = list()
        name_list = list()
        category_list = list()

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
            # columns = diff_pd.columns
            # fc_ind = list(columns).index('fc')
            # need_cols = ['seq_id', 'fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
            # need_cols += [columns[fc_ind - 2], columns[fc_ind - 1]]
            # print(need_cols)
            # samples = list()
            # for x in columns:
            #     _m = re.match(r'(.*)_count$', x)
            #     if _m:
            #         samples.append(_m.groups()[0])
            fname = os.path.basename(each)
            ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(diff_method.lower()), fname).groups()
            cmp_combine = ctrl + '|' + test
            name = '{}_vs_{}_volcano.pdf'.format(ctrl, test)
            field_index = str(pdf_name.index(name))
            cmp_list.append({'field':field_index, 'title':cmp_combine})

            # cmp_detail_dict[cmp_combine] = samples
            cmp_pd = pd.DataFrame([cmp_combine] * diff_pd.shape[0], columns=['compare'])
            # tmp_pd = pd.concat([diff_pd.loc[:, need_cols], cmp_pd], axis=1)
            # tmp_pd.columns = list(tmp_pd.columns[:-3]) + ['group1', 'group2', 'compare']
            # print(tmp_pd.columns)
            # diff_dict_list += tmp_pd.to_dict('records')

            # get volcano data
            status_list, stat_cutoff, down_num, up_num = self._get_volcano_status_cutoff(diff_pd, pvalue_padjust)
            # # sig_status[cmp_combine] = status_list
            # volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', pvalue_padjust, 'significant', 'regulate']]
            # bool_ind = volcano_pd[pvalue_padjust] <= 0
            # min_pvalue = min([x if x > 0 else '' for x in volcano_pd[pvalue_padjust].tolist()])
            # volcano_pd.loc[bool_ind, pvalue_padjust] = min_pvalue
            # volcano_pd[pvalue_padjust] = -np.log10(volcano_pd[pvalue_padjust])
            # volcano_pd.dropna(inplace=True)
            # regulate = list()
            # for i in volcano_pd.index:
            #     if volcano_pd['significant'].iloc[i] == 'no':
            #         regulate.append(status_list[0])
            #     if volcano_pd['significant'].iloc[i] == 'yes':
            #         if volcano_pd['regulate'].iloc[i] == 'up':
            #             regulate.append(status_list[2])
            #         if volcano_pd['regulate'].iloc[i] == 'down':
            #             regulate.append(status_list[1])
            #     # else:
            #     #     regulate.append(volcano_pd['regulate'].iloc[i])
            # volcano_pd['regulate'] = regulate
            # volcano_pd.columns = ['seq_id', 'log2fc', 'log10pvalue', 'significant', 'regulate']
            # bool_ind = volcano_pd['log10pvalue'] > stat_cutoff
            # volcano_pd.loc[bool_ind, 'log10pvalue'] = stat_cutoff
            # volcano_pd = pd.concat([volcano_pd, cmp_pd], axis=1)
            # volcano_pd_nosig = volcano_pd[volcano_pd['significant'] == 'no']
            # # random select 10000 not sig diff genes for plotting
            # if volcano_pd_nosig.shape[0] > 8000:
            #     volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
            #     volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.5)
            #     if volcano_pd_nosig.shape[0] > 12000:
            #         volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
            #     if volcano_pd_nosig.shape[0] > 12000:
            #         volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
            #     volcano_pd = pd.concat([volcano_pd_sig, volcano_pd_nosig], axis=0)
            # volcano_dict_list += volcano_pd.to_dict('records')

            # get column data
            value_list.append(down_num)
            name_list.append(cmp_combine)
            category_list.append('down')
            value_list.append(up_num)
            name_list.append(cmp_combine)
            category_list.append('up')




            # # get scatter data
            # scatter_pd = tmp_pd.loc[:, ['seq_id', 'group1', 'group2', 'compare', 'significant', 'regulate']]
            # scatter_pd.set_index('seq_id', inplace=True)
            # scatter_pd = scatter_pd.loc[volcano_pd['seq_id'], :].reset_index()
            # scatter_pd['group1'] = (scatter_pd['group1'] + 1).apply(np.log10)
            # scatter_pd['group2'] = (scatter_pd['group2'] + 1).apply(np.log10)
            # scatter_dict_list += scatter_pd.to_dict('records')
            # # get significant diff gene set and add to database
            # if create_geneset:
            #     name = ctrl + '_vs_' + test + '_' + exp_level
            #     # name += '_' + quant_method + '_' + diff_method
            #     # time_now = datetime.datetime.now()
            #     # name += '_' + time_now.strftime("%Y%m%d_%H%M%S")
            #     # name += '_' + time_now.strftime("%H%M%S")
            #     sig_seqs = list(diff_pd['seq_id'][diff_pd['significant'] == 'yes'])
            #     sig_regulate = list(diff_pd['regulate'][diff_pd['significant'] == 'yes'])
            #     if len(sig_seqs) >= 1:
            #         geneset_main_info = dict(
            #             project_sn=project_sn,
            #             task_id=task_id,
            #             version="v3",
            #             name=name,
            #             type=exp_level,
            #             desc='differential expressed gene set',
            #             group_id=group_id,
            #             source="diff_exp",
            #             gene_length=len(sig_seqs),
            #             is_use=1
            #         )
            #
            #         genet_detail_info = [{"seq_list": sig_seqs, "regulate_list": sig_regulate}]
            #         self.add_set(geneset_main_info, genet_detail_info)
        else:
            'loop end of diff_files'

        if main_id is None:
            # add main table info
            name = "DiffExpress" + '_' + quant_method + '_' + diff_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='differential analysis main table',
                # exp_level=exp_level,
                # group_dict=group_dict,
                # group_id=group_id,
                # control_id=cmp_id,
                # exp_id=exp_id,
                # diff_sum=dict(zip(levels[0][labels[0]], levels[1][labels[1]])),
                params=params,
                status="start",
            )
            main_id = self.create_db_table('{}_diff'.format(diff_method), [main_info])
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        # update main table
        # scatter_data = dict(name="seq_id",data=["log2fc", "log10pvalue"], category="regulate", condition={"type": "scatter"})
        # scatter_data = json.dumps(scatter_data)
        column_data = dict(name='name', data='value', category='category', condition={"type": "column"})
        column_data = json.dumps(column_data)
        pic_data = s3_output
        cmp_list = json.dumps(cmp_list)
        self.update_db_record('{}_diff'.format(diff_method), main_id,
                              cmp_list=cmp_list,
                              # cmp_detail=cmp_detail_dict,
                              # sig_status=sig_status,
                              # scatter_data=scatter_data,
                              column_data=column_data,
                              pic_data=pic_data,
                              main_id=main_id
                              )
        # # get id2gene_name dict
        # id2gene_name, id2desc, id2gid = self.get_gene_name_dict(main_id, 'sg_diff')
        #
        # # add detail info
        # for each in diff_dict_list:
        #     # each['gene_name'] = id2gene_name[each['seq_id']]
        #     # each['description'] = id2desc[each['seq_id']]
        #     if id2gid:
        #         each['gene_id'] = id2gid[each['seq_id']]
        # self.create_db_table('{}_diff_detail'.format(diff_method), diff_dict_list, tag_dict={'diff_id': main_id})

        # add summary detail
        # with open(summary, "r") as f:
        #     lines = f.readlines()
        #     if len(lines) >= 3:
        #         summary_pd = pd.read_table(summary, header=[0, 1])
        #         levels = summary_pd.columns.levels
        #         labels = summary_pd.columns.labels
        #         summary_pd.columns = levels[0][labels[0]]
        #         # summary_pd['gene_name'] = [id2gene_name[x] for x in summary_pd['seq_id']]
        #         # summary_pd['description'] = [id2desc[x] for x in summary_pd['seq_id']]
        #         # if id2gid:
        #         #     summary_pd['gene_id'] = [id2gid[x] for x in summary_pd['seq_id']]
        #         summary_dict_list = summary_pd.to_dict('records')
        #         self.create_db_table('{}_diff_summary'.format(diff_method), summary_dict_list, tag_dict={'diff_id': main_id})

                # add volcano detail
                # for each in volcano_dict_list:
                #     each['gene_name'] = id2gene_name[each['seq_id']]
        # self.create_db_table('{}_diff_volcano'.format(diff_method), volcano_dict_list, tag_dict={'diff_id': main_id})

                # # add scatter detail
                # for each in scatter_dict_list:
                #     each['gene_name'] = id2gene_name[each['seq_id']]
                # self.create_db_table('diff_scatter', scatter_dict_list, tag_dict={'diff_id': main_id})
            # else:
            #     pass

        column_df = pd.DataFrame({'value':value_list, 'name':name_list, 'category':category_list})
        column_df['type'] = 'column'
        column_df['diff_id'] = main_id
        self.create_db_table('{}_diff_column'.format(diff_method), column_df.to_dict('r'))
        # update status
        self.update_db_record('{}_diff'.format(diff_method), main_id, status="end", main_id=main_id, )
        return main_id

    def add_diffexp_noiseq(self, diff_output, exp_id=None, group_dict=None, group_id=None,
                    exp_level='T', quant_method='RSEM', diff_method='DESeq2', main_id=None,
                    project_sn='denovo_rna_v2', task_id='denovo_rna_v2', params=None,
                           create_geneset=True, nosig_info=None, s3_output=None, pdf_name=None):
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
            # self.create_db_table('sg_diff_detail', nosig_list, tag_dict={'diff_id': main_id, "seq_type": "no_sig"})
        summary = glob.glob(os.path.join(diff_output, '*diff_summary*.xls'))[0]
        # diff_method = os.path.basename(summary).split('_diff_summary.xls')[0]
        diff_files = glob.glob(os.path.join(diff_output, '*_vs_*.*.xls'))
        if not diff_files:
            self.bind_object.set_error('No target file found in %s', variables=(diff_output), code="53700204")
        cmp_list = list()
        cmp_detail_dict = dict()
        diff_dict_list = list()
        volcano_dict_list = list()
        sig_status = dict()
        scatter_dict_list = list()
        value_list = list()
        name_list = list()
        category_list = list()
        for each in diff_files:
            if each.endswith('.annot.xls'):
                continue
            if each.endswith('.normalize.xls'):
                continue
            if each.endswith('.sizeFactor.xls'):
                continue
            if each.endswith('.normFactor.xls'):
                continue
            fname = os.path.basename(each)
            ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(diff_method.lower()), fname).groups()
            diff_pd = pd.read_table(each, header=0, sep='\t')
            # columns = diff_pd.columns
            # fc_ind = list(columns).index('fc')
            # need_cols = ['seq_id', '{}_mean'.format(ctrl), '{}_mean'.format(test), 'fc', 'log2fc', 'D', 'prob',
            #              'significant', 'regulate']
            # need_cols += [columns[fc_ind - 4], columns[fc_ind - 3]]
            # print(need_cols)
            # samples = list()
            # for x in columns:
            #     _m = re.match(r'(.*)_count$', x)
            #     if _m:
            #         samples.append(_m.groups()[0])

            cmp_combine = ctrl + '|' + test
            # cmp_list.append(cmp_combine)
            name = '{}_vs_{}_volcano.pdf'.format(ctrl, test)
            field_index = str(pdf_name.index(name))
            cmp_list.append({'field':field_index, 'title':cmp_combine})
            # cmp_detail_dict[cmp_combine] = samples
            cmp_pd = pd.DataFrame([cmp_combine] * diff_pd.shape[0], columns=['compare'])
            # tmp_pd = pd.concat([diff_pd.loc[:, need_cols], cmp_pd], axis=1)
            # tmp_pd.columns = list(tmp_pd.columns[:-3]) + ['group1', 'group2', 'compare']
            # print(tmp_pd.columns)
            # diff_dict_list += tmp_pd.to_dict('records')

            # get volcano data
            status_list, down_num, up_num = self._get_volcano_status_cutoff_noiseq(diff_pd)
            # sig_status[cmp_combine] = status_list
            # volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', 'D', 'significant', 'regulate']]
            # volcano_pd.dropna(inplace=True)
            # regulate = list()
            # for i in volcano_pd.index:
            #     if volcano_pd['significant'].iloc[i] == 'no':
            #         regulate.append(status_list[0])
            #     if volcano_pd['significant'].iloc[i] == 'yes':
            #         if volcano_pd['regulate'].iloc[i] == 'up':
            #             regulate.append(status_list[2])
            #         if volcano_pd['regulate'].iloc[i] == 'down':
            #             regulate.append(status_list[1])
            # volcano_pd['regulate'] = regulate
            # volcano_pd.columns = ['seq_id', 'log2fc', 'D', 'significant', 'regulate']
            # volcano_pd = pd.concat([volcano_pd, cmp_pd], axis=1)
            # volcano_pd_nosig = volcano_pd[volcano_pd['significant'] == 'no']
            # # random select 10000 not sig diff genes for plotting
            # if volcano_pd_nosig.shape[0] > 8000:
            #     volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
            #     volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.5)
            #     if volcano_pd_nosig.shape[0] > 12000:
            #         volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
            #     if volcano_pd_nosig.shape[0] > 12000:
            #         volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
            #     volcano_pd = pd.concat([volcano_pd_sig, volcano_pd_nosig], axis=0)
            # volcano_dict_list += volcano_pd.to_dict('records')

            # get column data
            value_list.append(down_num)
            name_list.append(cmp_combine)
            category_list.append('down')
            value_list.append(up_num)
            name_list.append(cmp_combine)
            category_list.append('up')





            # # get scatter data
            # scatter_pd = tmp_pd.loc[:, ['seq_id', 'compare', 'significant', 'regulate', '{}_mean'.format(ctrl), '{}_mean'.format(test)]]
            # scatter_pd.set_index('seq_id', inplace=True)
            # scatter_pd = scatter_pd.loc[volcano_pd['seq_id'], :].reset_index()
            # scatter_pd['group1'] = (scatter_pd['{}_mean'.format(ctrl)]+1).apply(np.log10)
            # scatter_pd['group2'] = (scatter_pd['{}_mean'.format(test)]+1).apply(np.log10)
            # scatter_dict_list += scatter_pd.to_dict('records')
            # # get significant diff gene set and add to database
            # if create_geneset:
            #     name = ctrl + '_vs_' + test + '_' + exp_level
            #     # name += '_' + quant_method + '_' + diff_method
            #     # time_now = datetime.datetime.now()
            #     # name += '_' + time_now.strftime("%Y%m%d_%H%M%S")
            #     # name += '_' + time_now.strftime("%H%M%S")
            #     sig_seqs = list(diff_pd['seq_id'][diff_pd['significant'] == 'yes'])
            #     sig_regulate = list(diff_pd['regulate'][diff_pd['significant'] == 'yes'])
            #     if len(sig_seqs) >= 1:
            #         geneset_main_info = dict(
            #             project_sn=project_sn,
            #             task_id=task_id,
            #             version="v3",
            #             name=name,
            #             type=exp_level,
            #             desc='differential expressed gene set',
            #             group_id=group_id,
            #             source="diff_exp",
            #             gene_length=len(sig_seqs),
            #             is_use=1
            #         )
            #
            #         genet_detail_info = [{"seq_list": sig_seqs, "regulate_list": sig_regulate}]
            #         self.add_set(geneset_main_info, genet_detail_info)
        else:
            'loop end of diff_files'

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
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='differential analysis main table',
                exp_level=exp_level,
                group_dict=group_dict,
                group_id=group_id,
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
        # scatter_data = dict(name="seq_id",data=["log2fc", "D"], category="regulate", condition={"type": "scatter"})
        # scatter_data = json.dumps(scatter_data)
        column_data = dict(name='name', data='value', category='category', condition={"type": "column"})
        column_data = json.dumps(column_data)
        pic_data = s3_output
        cmp_list = json.dumps(cmp_list)
        self.update_db_record('{}_diff'.format(diff_method), main_id,
                              cmp_list=cmp_list,
                              # cmp_detail=cmp_detail_dict,
                              # sig_status=sig_status,
                              # scatter_data=scatter_data,
                              column_data=column_data,
                              pic_data=pic_data,
                              main_id=main_id
                              )
        # # get id2gene_name dict
        # id2gene_name, id2desc, id2gid = self.get_gene_name_dict(main_id, 'sg_diff')
        #
        # # add detail info
        # for each in diff_dict_list:
        #     # each['gene_name'] = id2gene_name[each['seq_id']]
        #     # each['description'] = id2desc[each['seq_id']]
        #     if id2gid:
        #         each['gene_id'] = id2gid[each['seq_id']]
        # self.create_db_table('sg_diff_detail', diff_dict_list, tag_dict={'diff_id': main_id})
        #
        # # add summary detail
        # with open(summary, "r") as f:
        #     lines = f.readlines()
        #     if len(lines) >= 3:
        #         summary_pd = pd.read_table(summary, header=[0, 1])
        #         levels = summary_pd.columns.levels
        #         labels = summary_pd.columns.labels
        #         summary_pd.columns = levels[0][labels[0]]
        #         # summary_pd['gene_name'] = [id2gene_name[x] for x in summary_pd['seq_id']]
        #         # summary_pd['description'] = [id2desc[x] for x in summary_pd['seq_id']]
        #         if id2gid:
        #             summary_pd['gene_id'] = [id2gid[x] for x in summary_pd['seq_id']]
        #         summary_dict_list = summary_pd.to_dict('records')
        #         self.create_db_table('sg_diff_summary', summary_dict_list, tag_dict={'diff_id': main_id})
        #
        #         # add volcano detail
        #         for each in volcano_dict_list:
        #             each['gene_name'] = id2gene_name[each['seq_id']]
        # self.create_db_table('{}_diff_volcano'.format(diff_method), volcano_dict_list, tag_dict={'diff_id': main_id})

        column_df = pd.DataFrame({'value':value_list, 'name':name_list, 'category':category_list})
        column_df['type'] = 'column'
        column_df['diff_id'] = main_id
        self.create_db_table('{}_diff_column'.format(diff_method), column_df.to_dict('r'))
            #     # add scatter detail
            #     for each in scatter_dict_list:
            #         each['gene_name'] = id2gene_name[each['seq_id']]
            #     self.create_db_table('sg_diff_scatter', scatter_dict_list, tag_dict={'diff_id': main_id})
            # else:
            #     pass

        # update status
        self.update_db_record('sg_diff', main_id, status="end", main_id=main_id, )
        return main_id
