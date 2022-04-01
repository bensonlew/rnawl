# -*- coding: utf-8 -*-
# __author__ = 'litangjian'
from collections import OrderedDict
import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import pandas as pd
import json
from mbio.api.database.denovo_rna_v2.api_base import ApiBase

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class TfApi(ApiBase):
    def __init__(self, bind_object):
        super(TfApi, self).__init__(bind_object)
        self._project_type = 'denovo_rna_v2'

    @report_check
    def add_tf_unigene(self, tf_unigene_path, bedpath=None, name=None, params=None, project_sn='denovo_rna_v2', task_id='denovo_rna_v2'):
        # bedpath是orf这个tool的output_dir里面的bed文件
        # task_id = self.bind_object.sheet.id
        # project_sn = self.bind_object.sheet.project_sn
        # 这个task_id，task_id也是自己随便写
        # task_id = "denovo_rna_v2"
        # project_sn = "denovo_rna_v2"
        tf_unigene = pd.read_table(tf_unigene_path, sep = "\t", header = 0)['Family'].value_counts().to_dict()
        # params = {"search_pfam": 'True', "p_length":50, "Markov_length": 3000, "E": 1e-3, "cpu": 20, "hmmcan1": "noali", "hmmcan2": "acc", "hmmcan3": "notextw"}
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'TF_' + "G" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'version': "v2",
            'status': 'start',
            'desc': 'tf_unigene结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'tf':tf_unigene,
            'type':'G',
        }
        tf_unigene_id = self.create_db_table('sg_tf', [insert_data])
        self. add_tf_unigene_detail(tf_unigene_path, tf_unigene_id)
        self.update_db_record('sg_tf', tf_unigene_id, status="end", main_id=tf_unigene_id)
        self.update_sgtask_record(task_id=task_id, bedpath=bedpath)
        return tf_unigene_id

    @report_check
    def add_tf_unigene_detail(self, tf_unigene_path, tf_unigene_id, name=None):
        # task_id = self.bind_object.sheet.id
        # project_sn = self.bind_object.sheet.project_sn
        # 这个task_id，task_id也是自己随便写
        # task_id = "denovo_rna_v2"
        # project_sn = "denovo_rna_v2"
        tf_unigene = pd.read_table(tf_unigene_path, sep = "\t", header = 0)['Family'].value_counts().to_dict()
        data_list = []
        # params = {"search_pfam": 'True', "p_length":50, "Markov_length": 3000, "E": 1e-3, "cpu": 20, "hmmcan1": "noali", "hmmcan2": "acc", "hmmcan3": "notextw"}
        # params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        with open(tf_unigene_path, 'r') as f1:
            f1.readline()
            for line in f1:
                line = line.strip().split('\t')
                data = [
                    ('tf_id', tf_unigene_id),
                    ('target_name', line[0]),
                    ('accession', line[1]),
                    ('evalue', float(line[3])),
                    ('score', float(line[4])),
                    ('keyword', line[6] + '|' + line[7] + '|' + line[5]),
                    ('seq_id', line[8]),
                    # ('type', 'G')
                    ]
                data = SON(data)
                data_list.append(data)

        self.create_db_table('sg_tf_detail', data_list)

        data_list = []
        with open(tf_unigene_path + ".stat", 'r') as f1:
            f1.readline()
            for line in f1:
                line = line.strip().split('\t')
                data = [
                    ('tf_id', tf_unigene_id),
                    ('tf_family', line[0]),
                    ('num', line[1]),
                    ('seq_list', line[2].split(","))
                    ]
                data = SON(data)
                data_list.append(data)

        self.create_db_table('sg_tf_stat', data_list)


    @report_check
    def add_tf_transcript(self, tf_transcript_path, bedpath=None, name=None, params=None, project_sn='denovo_rna_v2', task_id='denovo_rna_v2'):
        # task_id = self.bind_object.sheet.id
        # project_sn = self.bind_object.sheet.project_sn
        # 这个task_id，task_id也是自己随便写
        # task_id = "denovo_rna_v2"
        # project_sn = "denovo_rna_v2"
        tf_transcript = pd.read_table(tf_transcript_path, sep="\t", header=0)['Family'].value_counts().to_dict()
        # params = {"search_pfam": 'True', "p_length":50, "Markov_length": 3000, "E": 1e-3, "cpu": 20, "hmmcan1": "noali", "hmmcan2": "acc", "hmmcan3": "notextw"}
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'TF_' + "T"+ str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'version': "v2",
            'status': 'start',
            'desc': 'tf_transcript结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'tf': tf_transcript,
            'type':'T'

        }
        tf_transcript_id = self.create_db_table('sg_tf', [insert_data])
        self.add_tf_transcript_detail(tf_transcript_path, tf_transcript_id)
        self.update_db_record('sg_tf', tf_transcript_id, status="end", main_id=tf_transcript_id)
        # self.update_sgtask_record(task_id=task_id, bedpath=bedpath)
        return tf_transcript_id

    def add_tf_transcript_detail(self, tf_transcript_path, tf_transcript_id, name=None):
        # task_id = self.bind_object.sheet.id
        # project_sn = self.bind_object.sheet.project_sn
        # 这个task_id，task_id也是自己随便写
        # task_id = "denovo_rna_v2"
        # project_sn = "denovo_rna_v2"
        tf_transcript = pd.read_table(tf_transcript_path, sep="\t", header=0)['Family'].value_counts().to_dict()
        data_list = []
        with open(tf_transcript_path, 'r') as f1:
            f1.readline()
            for line in f1:
                line = line.strip().split('\t')
                data = [
                    ('tf_id', tf_transcript_id),
                    ('target_name', line[0]),
                    ('accession', line[1]),
                    ('seq_id', line[2]),
                    ('evalue', float(line[3])),
                    ('score', float(line[4])),
                    ('keyword', line[6] + '|' + line[7] + '|' + line[5]),
                    # ('type', 'T')
                ]
                data = SON(data)
                data_list.append(data)

        self.create_db_table('sg_tf_detail', data_list)

        data_list = []

        with open(tf_transcript_path  + ".stat", 'r') as f1:
            f1.readline()
            for line in f1:
                line = line.strip().split('\t')
                data = [
                    ('tf_id', tf_transcript_id),
                    ('tf_family', line[0]),
                    ('num', line[1]),
                    ('seq_list', line[2].split(","))
                    ]
                data = SON(data)
                data_list.append(data)

        self.create_db_table('sg_tf_stat', data_list)


if __name__ == '__main__':
    anno = TfApi(None)

    from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
    from biocluster.wsheet import Sheet
    import random


    data = {
        "id": "denovo_rna_v2_upgrade",
        #+ str(random.randint(1,10000)),
        #"id": "denovo_rna_v2",
        "project_sn": "denovo_rna_v2_upgrade",
        #+ str(random.randint(1,10000)),
        "type": "workflow",
        "name": "denovo_rna_v2.denovo_test_api",
        "options": {
        },
    }
    wsheet = Sheet(data=data)
    wf = DenovoTestApiWorkflow(wsheet)

    wf.IMPORT_REPORT_DATA = True
    wf.IMPORT_REPORT_AFTER_END = False
    wf.test_api = wf.api.api("denovo_rna_v2.tf_api")

    tf_unigene_path = '/mnt/ilustre/users/sanger-dev/workspace/20190807/Single_annot_orf15-19-50/AnnotOrfpfam/Predict/merge_only_unigene_animal'
    bedpath = '/mnt/ilustre/users/sanger-dev/workspace/20190807/Single_annot_orf15-19-50/AnnotOrfpfam/Predict/all_predicted.bed'
    wf.test_api.add_tf_unigene(tf_unigene_path, bedpath=bedpath, name=None, params=None, project_sn='denovo_rna_v2_upgrade', task_id='denovo_rna_v2_upgrade')
    tf_transcript_path = '/mnt/ilustre/users/sanger-dev/workspace/20190807/Single_annot_orf15-19-50/AnnotOrfpfam/Predict/merge_only_transcript_animal'
    wf.test_api.add_tf_transcript(tf_transcript_path, name=None, params=None, project_sn='denovo_rna_v2_upgrade', task_id='denovo_rna_v2_upgrade')
