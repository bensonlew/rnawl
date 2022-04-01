# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import datetime
import json
import os
import types
import unittest

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from api_base import ApiBase


class SmallRnaMapping(ApiBase):  # 继承父类Apibase
    def __init__(self, bind_object):
        super(SmallRnaMapping, self).__init__(bind_object)
        self.task_id = self.bind_object.sheet.id
        self.rock_params = json.dumps({"software": "Rockhopper"}, sort_keys=True, separators=(',', ':'))
        self.bind_object.logger.info("********rock_param is: {} ***************".format(self.rock_params))
        # self.promote_params = json.dumps({""}, sort_keys=True, separators=(',', ':'))

    @report_check  # 比对分析结果---主表
    def add_mapping(self, sample_list=None, result_dir=None, name=None, params=None):
        """
        主表
        smallrna 比对结果统计表
        """
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        params['samples'] = ",".join(sample_list)  ### 将列表变为字符串
        self.bind_object.logger.info("#############  params['samples']: {} ############### ".format(params['samples']))
        self.bind_object.logger.info("#############  params is: {}".format(
            params))  # params is: {'chr_num': '10', 'task_type': '2', 'samples': 'A1,A2,C2,C3,S1,S3'} 是个字典。

        if params is None:
            params_dict = dict()
        elif type(params) == dict:
            params_dict = params
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))  # 将字典转化为字符串
        # params= "{\"chr_num\":\"10\",\"samples\":\"A1,A2,C2,C3,S1,S3\",\"task_type\":\"2\"}",
        else:
            params_dict = json.loads(params)

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'SmallRnaMapping_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            # 'result_dir': result_dir,
            'status': 'start',  # 用于判断数据已经插入mongo库中，"end" 表示已成功插入！
            'desc': '比对分析结果主表',
            'type': "origin",
            'sample_list': sample_list,
            'sample_list_circos': sample_list,
            'result_dir': result_dir,
            'result_origin': result_dir,
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }

        collection = self.db['mapping']
        mapping_id = collection.insert_one(
            insert_data).inserted_id  # collection.insert_one(insert_data) 产生一个对象, .inserted_id 返回对象的属性。
        self.bind_object.logger.info("已在mongo库中创建collection：mapping")
        return mapping_id

    @report_check
    def add_circos(self, sample_list=None, result_dir=None, name=None, params=None):  # 主表
        """
        主表
        smallrna 圈图数据
        """
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        params['samples'] = ",".join(sample_list)

        if params is None:
            params_dict = dict()
        elif type(params) == dict:
            params_dict = params
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))  # 将字典转化为字符串
        # params= "{\"chr_num\":\"10\",\"samples\":\"A1,A2,C2,C3,S1,S3\",\"task_type\":\"2\"}",
        else:
            params_dict = json.loads(params)

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'SmallRnaCircos_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            # 'result_dir': result_dir,
            'status': 'start',
            'desc': '圈图分析结果主表',
            'type': "origin",
            'sample_list': sample_list,
            'sample_list_circos': sample_list,
            'result_dir': result_dir,
            'result_origin': result_dir,
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'version': 'v1'
        }
        collection = self.db['circos']  # 创建collection，如果创建新的collection，那么这个mapping_id 也应相应改变。
        circos_id = collection.insert_one(insert_data).inserted_id  # 集合中插入文档使用 insert_one() 方法
        self.bind_object.logger.info("circos_id is:{}".format(circos_id))  # 输出circos id
        self.update_db_record(table_name='circos', record_id=circos_id, status="end", main_id=circos_id)  # 对表格记录进行更新
        self.bind_object.logger.info("已在mongo库中创建collection：circos")

    @report_check
    def add_map_stat(self, sample_list=None, map_dir=None, sample_file=None, params=None):
        """
        导入 rockhopper 得到的opera, utr tss 等信息
        fq_list = sample_file:      mirna_qc Module module 质控结果 output_dir文件夹
        map_dir = mapping_dir:      mapper_and_stat  Tool  比对结果 output_dir文件夹
        """
        if os.path.exists(map_dir):
            pass
        else:
            raise Exception('与基因组文件{} 不存在，请检查'.format(map_dir))

        # 以下是为了获得样本名
        if os.path.exists(sample_file):  # 也就是下面的 fq_list
            sample_list = list()
            with open(sample_file, 'r') as f:
                for line in f.readlines():
                    sample = line.strip().split("\t")[1]
                    if sample in sample_list:
                        pass
                    else:
                        sample_list.append(sample)
        self.bind_object.logger.info(
            "本次项目的样本名为：{}".format(sample_list))  # 本次项目的样本名为：['A1', 'A2', 'C2', 'C3', 'S1', 'S3']
        # 以上的遍历获取的是本项目的样本名，
        task_id = self.task_id
        self.remove_table_by_main_record(main_table='mapping', task_id=task_id, detail_table=['mapping_detail'],
                                         detail_table_key='mapping_id')
        mapping_id = self.add_mapping(sample_list, map_dir, params=params)  #
        # 需调用add_mapping方法，而该方法会创建两个collection，名叫"mapping"、"circos"。
        self.add_mapping_detail(mapping_id=mapping_id, map_dir=map_dir, sample_list=sample_list)
        # 需调用add_mapping_detail方法，该方法会创建一个collection，为"mapping_detail"，"mapping_stat"。
        self.update_db_record('mapping', mapping_id, status="end", main_id=mapping_id)

    @report_check  # 比对分析详情表（作图）
    def add_mapping_detail(self, mapping_id=None, map_dir=None, sample_list=None):
        if not isinstance(mapping_id, ObjectId):
            if isinstance(mapping_id, types.StringTypes):
                mappping_id = ObjectId(mapping_id)
            else:
                raise Exception('mapping_id必须为ObjectId对象或其对应的字符串！')
        data_stat_list = []
        data_list = []

        all_stat = pd.read_table(map_dir + "/All_map_stat.xls", header=0, index_col=0,
                                 usecols=['Chromosome', 'Total_num'])

        with open(map_dir + "/Genome_map_stat.xls", 'w') as f:
            f.write("Sample\tTotal reads\tMapped reads\tMapped reads(+)\tMapped reads(-)\n")

            all_stat = all_stat[2:]
            all_stat['ref'] = all_stat.index
            all_stat.rename(str.lower, axis='columns', inplace=True)
            for sample in sample_list:
                sample_stat = pd.read_table(map_dir + "/" + sample + "_map_stat.xls", header=0, index_col=0,
                                            usecols=['Chromosome', 'Total_num', 'Forword.1', 'Reverse.1'])[2:]
                all_stat[sample + '_total'] = sample_stat['Total_num'].map(int)
                all_stat[sample + '_for'] = sample_stat['Forword.1'].map(int)
                all_stat[sample + '_rev'] = sample_stat['Reverse.1'].map(int)

                stat2 = pd.read_table(map_dir + "/" + sample + "_map_stat.xls", header=0, index_col=0)[:2]
                data = dict({
                    "sample_name": sample,
                    "total_reads": int(stat2.iloc[0]['Total_num']),
                    "mapped_reads": int(stat2.iloc[1]['Total_num']),
                    "mapped_reads_for": int(stat2.iloc[1]['Forword.1']),
                    "mapped_reads_rec": int(stat2.iloc[1]['Reverse.1']),
                    "mapping_id": mapping_id
                })
                data_stat_list.append(data)

                f.write("{}\t{}\t{}\t{}\t{}\n".format(sample,
                                                      int(stat2.iloc[0]['Total_num']),
                                                      int(stat2.iloc[1]['Total_num']),
                                                      stat2.iloc[1]['Forword.1'],
                                                      stat2.iloc[1]['Reverse.1']))

        all_stat['mapping_id'] = mapping_id
        all_stat = all_stat.fillna(0)

        def float2int(x):
            if type(x) == "float":
                return int(x)
            else:
                return x

        all_stat = all_stat.applymap(float2int)

        columns = []
        columns.extend([x + "_total" for x in sample_list])

        columns_rename = []
        columns_rename.extend(sample_list)
        all_stat_choose = all_stat.loc[:, columns]
        rename_dict = dict(zip(columns, columns_rename))
        all_stat_choose = all_stat_choose.rename(columns=rename_dict)
        all_stat_choose.to_csv(map_dir + "/Chro_map_stat.xls", sep="\t", index=True)

        row_dict_list = all_stat.to_dict('records')
        main_collection = self.db['mapping_detail']

        try:
            self.create_db_table('mapping_detail', row_dict_list)
            self.create_db_table('mapping_stat', data_stat_list)
        except Exception as e:
            raise Exception("导入main: %s出错!" % (all_stat))


class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """

    def test_mongo(test):
        from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "small_rna_0929",
            "project_sn": "small_rna",
            "type": "workflow",
            "name": "whole_transcriptome.small_rna_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = SmallRnaTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.api_map = wf.api.api("whole_transcriptome.small_rna_mapping")
        # mapping_dir = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/data5/mapping_out/"
        mapping_dir = "/mnt/ilustre/users/sanger-dev/workspace/20190928/Single_MapperAndStat13-39-10/MapperAndStat/output"
        # fq_list = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/list.txt"
        fq_list = "/mnt/ilustre/users/sanger-dev/workspace/20190926/Single_whole_transcriptome.smallrna.mirna_qc15-22-40/MirnaQc/output/clean_data/list.txt"
        params = {
            "chr_num": "10",
            "task_type": "2"}

        wf.api_map.add_map_stat(sample_list=None, map_dir=mapping_dir, sample_file=fq_list, params=params)
        # wf.api_map.add_mapping_detail(mapping_id=mapping_id, map_dir=mapping_dir) # sample_list=None


if __name__ == '__main__':
    unittest.main()
