# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from collections import OrderedDict
import os
import json
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import pandas as pd
from collections import OrderedDict
from mbio.api.database.denovo_rna_v2.api_base import ApiBase
import numpy as np

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class SnpApi(ApiBase):
    def __init__(self, bind_object):
        super(SnpApi, self).__init__(bind_object)

    @report_check
    def add_snp(self,new_snp_rewrite=None, new_indel_rewrite=None, depth_path=None, hh_path=None, tt_new_per_path=None,
        call_vcf_path=None, name=None, params=None, project_sn='denovo_rna_v2', task_id='denovo_rna_v2', is_report='no'):
        # task_id = self.bind_object.sheet.id
        # project_sn = self.bind_object.sheet.project_sn
        # 这个task_id，task_id也是自己随便写
        # task_id = "denovo_rna_v2"
        # project_sn = "denovo_rna_v2"
        # df = pd.read_table(new_snp_rewrite, sep = "\t", header = 0, encoding = 'utf8')
        # indel_col_len = df.shape[1]
        # select_col = range(5, indel_col_len, 4)
        # df_selected = df.iloc[:,select_col].replace([r'-'], -1)
        # sample_names = list(df_selected.columns)
        # call_vcf_path = call_vcf_path
        # params = {"qual": 20, "dp":1}
        # params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        # insert_data = {
        #     'project_sn': project_sn,
        #     'task_id': task_id,
        #     'name': name if name else 'snp_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
        #     'params': params,
        #     'status': 'start',
        #     'desc': 'snp_结果表',
        #     'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        #     'sample_names': sample_names,
        #     'call_vcf_path': call_vcf_path,
        #     'is_report':is_report
        #            }
        #
        # snp_id = self.create_db_table('sg_snp', [insert_data])
        # return snp_id

        if new_snp_rewrite is None and new_indel_rewrite is None:
            # return snp_id
            self.logger.info("此次分析没有call出snp和indel")

        if new_snp_rewrite is not None and new_indel_rewrite is not None:
            df = pd.read_table(new_snp_rewrite, sep="\t", header=0, encoding='utf8')
            indel_col_len = df.shape[1]
            select_col = range(5, indel_col_len, 4)
            df_selected = df.iloc[:, select_col].replace([r'-'], -1)
            sample_names = list(df_selected.columns)
            call_vcf_path = call_vcf_path
            # params = {"qual": 20, "dp": 1}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'name': name if name else 'snp_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                'params': params,
                'status': 'start',
                'desc': 'snp_结果表',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'sample_names': sample_names,
                'call_vcf_path': call_vcf_path,
                'is_report': is_report
            }
            snp_id = self.create_db_table('sg_snp', [insert_data])

            self.add_snp_detail(new_snp_rewrite, snp_id, name=None, params=None)
            self.add_indel_detail(new_indel_rewrite, snp_id, name=None, params=None)
            self.add_depth_type(depth_path, snp_id)
            self.add_hh_type(hh_path, snp_id)
            self.add_tt_new_per_type(tt_new_per_path, snp_id)
            self.update_db_record('sg_snp', snp_id, status="end", main_id=snp_id)
            # return snp_id

        if new_snp_rewrite is not None and new_indel_rewrite is None:
            df = pd.read_table(new_snp_rewrite, sep="\t", header=0, encoding='utf8')
            indel_col_len = df.shape[1]
            select_col = range(5, indel_col_len, 4)
            df_selected = df.iloc[:, select_col].replace([r'-'], -1)
            sample_names = list(df_selected.columns)
            call_vcf_path = call_vcf_path
            # params = {"qual": 20, "dp": 1}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'name': name if name else 'snp_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                'params': params,
                'status': 'start',
                'desc': 'snp_结果表',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'sample_names': sample_names,
                'call_vcf_path': call_vcf_path,
                'is_report': is_report
            }

            snp_id = self.create_db_table('sg_snp', [insert_data])

            self.add_snp_detail(new_snp_rewrite, snp_id, name=None, params=None)
            # self.add_indel_detail(new_indel_rewrite, snp_id, name=None, params=None)
            self.add_depth_type(depth_path, snp_id)
            self.add_hh_type(hh_path, snp_id)
            self.add_tt_new_per_type(tt_new_per_path, snp_id)
            self.update_db_record('sg_snp', snp_id, status="end", main_id=snp_id)
            # return snp_id

        if new_snp_rewrite is None and new_indel_rewrite is not None:
            df = pd.read_table(new_indel_rewrite, sep="\t", header=0, encoding='utf8')
            indel_col_len = df.shape[1]
            select_col = range(5, indel_col_len, 4)
            df_selected = df.iloc[:, select_col].replace([r'-'], -1)
            sample_names = list(df_selected.columns)
            call_vcf_path = call_vcf_path
            # params = {"qual": 20, "dp": 1}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'name': name if name else 'snp_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                'params': params,
                'status': 'start',
                'desc': 'snp_结果表',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'sample_names': sample_names,
                'call_vcf_path': call_vcf_path,
                'is_report': is_report
            }

            snp_id = self.create_db_table('sg_snp', [insert_data])
            self.update_db_record('sg_snp', snp_id, status="end",main_id=snp_id)
            # self.add_snp_detail(new_snp_rewrite, snp_id, name=None, params=None)
            self.add_indel_detail(new_indel_rewrite, snp_id, name=None, params=None)
            # return snp_id

    @report_check
    def add_snp_detail(self,new_snp_rewrite, snp_id, name=None, params=None):
        data_list = []
        if os.path.isfile(new_snp_rewrite):
            df = pd.read_table(new_snp_rewrite, sep ="\t", header = 0)
            snp_col_len = df.shape[1]
            select_col = range(6, snp_col_len, 4)
            for col in select_col:
                tmp_col = df.iloc[:, col]
                tmp_list = [int(x) if str(x).isdigit() else x for x in tmp_col]
                df.iloc[:, col] = tmp_list
            df.columns = ['seq_id'] + list(df.columns)[1:]
            # df.columns = [x.lower() for x in list(df.columns)]
            df['type'] = "SNP"
            df['snp_id'] = ObjectId(snp_id)
            data_list = df.to_dict('records')
            self.create_db_table('sg_snp_detail', data_list)

    @report_check
    def add_indel_detail(self,new_indel_rewrite, snp_id, name=None, params=None):
        def str2int(listtype):
            length = len(listtype)
            list1 = list()
            for x in listtype[6:length:4]:
                if x == '-':
                    x = x
                    list.append(x)
                else:
                    x = int(x)
                    list.append(x)

        data_list = []
        if os.path.isfile(new_indel_rewrite):
            df = pd.read_table(new_indel_rewrite, sep ="\t", header = 0)
            snp_col_len = df.shape[1]
            select_col = range(6, snp_col_len, 4)
            for col in select_col:
                tmp_col = df.iloc[:, col]
                tmp_list = [int(x) if str(x).isdigit() else x for x in tmp_col]
                df.iloc[:, col] = tmp_list
            df.columns = ['seq_id'] + list(df.columns)[1:]
            # df.columns = [x.lower() for x in list(df.columns)]
            df['type'] = "InDel"
            df['snp_id'] = ObjectId(snp_id)
            data_list = df.to_dict('records')
            self.create_db_table('sg_snp_detail', data_list)

    def add_depth_type(self, depth_path, snp_id):
        data_list = []
        if os.path.isfile(depth_path):
            # with open(depth_path, 'r') as f1:
            #     header = f1.readline()
            #     head_list = [x.lower() for x in header.strip().split("\t")] + ['snp_id'] + ['stat_type']
            #     for line in f1:
            #         line = [int(x) if x.isdigit() else x for x in line.strip().split('\t')] + [snp_id] + ['depth']
            #         data_dcit = OrderedDict(zip(head_list, line))
            #         data_list.append(data_dcit)
            df = pd.read_table(depth_path, sep ="\t", header = 0)
            df = pd.read_table(depth_path, sep ="\t", header = 0, )
            # df.columns = [x.lower() for x in list(df.columns)]
            df.set_index('depth', inplace=True)
            df.loc['total'] = df.sum(axis=0)
            df.reset_index(level=0, inplace=True)
            df['snp_id'] = ObjectId(snp_id)
            df.fillna("", inplace=True)
            # records会导致整数变为小数类型，就是保持数据的类型变为一致
            data_list = df.to_dict('records')
            for tmp_dict in data_list:
                for k in tmp_dict:
                    if k.endswith('_depth'):
                        tmp_dict[k] = int(tmp_dict[k])

            self.create_db_table('sg_snp_dp', data_list)

    def add_hh_type(self, hh_path, snp_id):
        data_list = []
        if os.path.isfile(hh_path):
            with open(hh_path, 'r') as f1:
                header = f1.readline()
                head_list = [x.lower() for x in header.strip().split("\t")] + ['snp_id']
                # head_list = [x.lower() for x in header.strip().split("\t")] + ['snp_id']
                for line in f1:
                    line = [int(x) if x.isdigit() else x for x in line.strip().split('\t')] + [ObjectId(snp_id)]
                    data_dcit = OrderedDict(zip(head_list, line))
                    data_list.append(data_dcit)

            self.create_db_table('sg_snp_hh', data_list)

    def add_tt_new_per_type(self, tt_new_per_path, snp_id):
        data_list = []
        if os.path.isfile(tt_new_per_path):
            df = pd.read_table(tt_new_per_path, sep ="\t", header = 0)
            # df.columns = [x.lower() for x in list(df.columns)]
            df['type'] = [x[0] + '/' + x[1] for x in list(df['type'])]
            df.set_index('type', inplace=True)
            df.loc['total'] = df.sum(axis=0)
            df.reset_index(level=0, inplace=True)
            df['snp_id'] = ObjectId(snp_id)
            df.fillna("", inplace=True)
            data_list = df.to_dict('records')
            for tmp_dict in data_list:
                for k in tmp_dict:
                    if not k.endswith('_per') and not k.endswith('_id') and k != 'type':
                        tmp_dict[k] = int(tmp_dict[k])
            self.create_db_table('sg_snp_tt', data_list)

if __name__ == '__main__':
    anno = SnpApi(None)
    new_snp_rewrite = '/mnt/ilustre/users/sanger-dev/workspace/20171228/Single_module_snp34235/Snp/Snpfinal/new_snp_rewrite'
    depth_path = '/mnt/ilustre/users/sanger-dev/workspace/20171228/Single_module_snp34235/Snp/Snpfinal/depth_new_per'
    hh_path = '/mnt/ilustre/users/sanger-dev/workspace/20171228/Single_module_snp34235/Snp/Snpfinal/statis_hh'
    tt_new_per_path = '/mnt/ilustre/users/sanger-dev/workspace/20171228/Single_module_snp34235/Snp/Snpfinal/tt_new_per'
    call_vcf_path = "/mnt/ilustre/users/sanger-dev/workspace/20171202/Single_test_snp40206/Snp/call.vcf"
    new_indel_rewrite = None
    anno.add_snp(new_snp_rewrite, new_indel_rewrite, depth_path, hh_path, tt_new_per_path, call_vcf_path, name=None, params=None)
    # print snp_id
    # anno.add_snp_detail(new_snp_rewrite, name=None, params=None)
    # anno.add_depth_type(depth_path)
    # anno.add_hh_type(hh_path)
    # anno.add_tt_new_per_type(tt_new_per_path)
