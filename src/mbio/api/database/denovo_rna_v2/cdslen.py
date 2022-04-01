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
import json
from api_base import ApiBase
import pandas as pd

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Cdslen(ApiBase):
    def __init__(self, bind_object):
        super(Cdslen, self).__init__(bind_object)
        self._project_type = 'denovo_rna_v2'

    @report_check
    def add_cds_unigene_length(self, cds_unigene_length, name=None, params=None, project_sn='denovo_rna_v2_upgrade', task_id='denovo_rna_v2_upgrade'):
        # params = {"search_pfam": 'True', "p_length":50, "Markov_length": 3000, "E": 1e-3, "cpu": 20, "hmmcan1": "noali", "hmmcan2": "acc", "hmmcan3": "notextw"}
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        cds_dcit = OrderedDict()
        with open(cds_unigene_length, 'r') as f2:
            _ = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                cds_dcit[line[0]] = int(line[1])

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'Cdslen_' + "G"+ str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'start',
            'desc': 'cds_unigene_length结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'cds_dcit':cds_dcit,
            'type':'G'
        }
        cds_unigene_id = self.create_db_table('sg_cds', [insert_data])
        # self.update_db_record('sg_cds', cds_unigene_id, status="end", main_id=cds_unigene_id)
        record_dict = {"_id": cds_unigene_id, "task_id": task_id}
        self.update_db_record_by_dict('sg_cds', record_dict, status="end", main_id=cds_unigene_id)
        return cds_unigene_id

    @report_check
    def add_cds_transcript_length(self, cds_transcript_length, name=None, params=None, project_sn='denovo_rna_v2_upgrade', task_id='denovo_rna_v2_upgrade'):
        # params = {"search_pfam": 'True', "p_length":50, "Markov_length": 3000, "E": 1e-3, "cpu": 20, "hmmcan1": "noali", "hmmcan2": "acc", "hmmcan3": "notextw"}
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        cds_dcit = OrderedDict()
        with open(cds_transcript_length, 'r') as f2:
            _ = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                cds_dcit[line[0]] = int(line[1])

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'Cdslen_' + "T" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'start',
            'desc': 'cds_transcript_length结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'cds_dcit': cds_dcit,
            'type':'T'
        }
        cds_transcript_id = self.create_db_table('sg_cds', [insert_data])
        record_dict = {"_id": cds_transcript_id, "task_id": task_id}
        self.update_db_record_by_dict('sg_cds', record_dict, status="end", main_id=cds_transcript_id)
        #self.update_db_record('sg_cds', cds_transcript_id, status="end", main_id=cds_transcript_id)
        return cds_transcript_id

    @report_check
    def cds(self, cds_detail, tran2unigene=None, name=None, params=None, project_sn='denovo_rna_v2_upgrade', task_id='denovo_rna_v2_upgrade', cds_unigene_length = None, cds_transcript_length =None):
        project_sn = self.bind_object.sheet.project_sn
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        cds_unigene_length_dict = OrderedDict()
        cds_transcript_length_dict = OrderedDict()
        with open(cds_unigene_length, 'r') as f2:
            _ = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                cds_unigene_length_dict[line[0]] = int(line[1])

        with open(cds_transcript_length, 'r') as f2:
            _ = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                cds_transcript_length_dict[line[0]] = int(line[1])


        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'TF_' + "G" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'version': "v2",
            "unigene_len": cds_unigene_length_dict,
            "transcript_len": cds_transcript_length_dict,
            'status': 'start',
            'desc': 'cds预测结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),

        }
        cds_id = self.create_db_table('sg_cds', [insert_data])
        self.add_cds_detail(cds_id, cds_detail, tran2unigene)
        # self.update_db_record('sg_cds', cds_id, status="end", main_id=cds_id)
        # 因为数据库分片修改
        self.update_db_record_by_dict('sg_cds', {"task_id": task_id}, status="end", main_id=cds_id)


    @report_check
    def add_cds_detail(self, cds_id, cds_detail, tran2unigene):
        cds_pd = pd.read_table(cds_detail, header=0)
        cds_pd['cds_id'] = ObjectId(cds_id)
        cds_pd['is_gene'] = cds_pd['is_gene'].map(lambda x:True if x == "yes" else False)
        result_dict_list = cds_pd.to_dict("records")
        self.create_db_table('sg_cds_detail', result_dict_list)


if __name__ == '__main__':
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
    wf.test_api = wf.api.api("denovo_rna_v2.cdslen")


    cds_unigene_length_path = '/mnt/ilustre/users/sanger-dev/workspace/20190807/Single_annot_orf15-19-50/AnnotOrfpfam/Predict/output/cds_len_unigene.txt'
    wf.test_api.add_cds_unigene_length(cds_unigene_length_path, name=None, params=None)
    cds_transcript_length_path = '/mnt/ilustre/users/sanger-dev/workspace/20190807/Single_annot_orf15-19-50/AnnotOrfpfam/Predict/output/cds_len_transcript.txt'
    wf.test_api.add_cds_transcript_length(cds_transcript_length_path, name=None, params=None)
    wf.test_api.cds("/mnt/ilustre/users/sanger-dev/workspace/20190807/Single_annot_orf15-19-50/AnnotOrfpfam/Predict/output/all_predicted.xls", task_id="denovo_rna_v2_upgrade", cds_unigene_length= cds_unigene_length_path, cds_transcript_length=  cds_transcript_length_path)
