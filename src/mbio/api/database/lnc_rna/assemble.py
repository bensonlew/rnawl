# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'

from mbio.api.database.lnc_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import os
import datetime
import re
import json
from bson.objectid import ObjectId
from bson.son import SON
import unittest

class Assemble(ApiBase):
    def __init__(self, bind_object):
        super(Assemble, self).__init__(bind_object)
        self._project_type = 'lnc_rna'

    @report_check
    def add_assemble_result(self, params=None, all_gtf_path=None, merged_path=None, statistics_path=None):
        self.bind_object.logger.info('start creating table in sg_transcripts')
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        merged_list = list()
        merged_gtf_path = dict()
        merged_gtf_path['gtf'] = os.path.join(merged_path, 'merged.gtf')
        merged_list.append(merged_gtf_path)
        merged_fa_path = dict()
        merged_fa_path['fa'] = os.path.join(merged_path, 'merged.fa')
        merged_list.append(merged_fa_path)
        time_now = datetime.datetime.now()
        name = 'Assemble_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'start',
            'desc': 'assemble_main_table',
            'created_ts': created_ts,
            'step': [200, 300, 600, 1000],
            'seq_type': ['i', 'j', 'o', 'u', 'x'],
            'step_of_trans_or_genes': [1, 5, 10, 20],
        }
        collection = self.db['sg_transcripts']
        transcript_id = self.create_db_table('sg_transcripts', [insert_data])
        self.bind_object.logger.info('succeed in creating table in sg_transcripts')
        code_files = os.path.join(statistics_path, 'code_num.txt')
        self.add_transcripts_step(transcript_id=transcript_id, statistics_path=statistics_path)
        self.add_transcripts_relations(transcript_id=transcript_id, statistics_path=statistics_path)
        self.add_transcripts_seq_type(transcript_id=transcript_id, code_file=code_files)
        self.update_db_record('sg_transcripts', transcript_id, status='end', main_id=transcript_id)
        self.bind_object.logger.info('succeed in updating table in sg_transcripts')
        return transcript_id

    def add_transcripts_step(self, transcript_id, statistics_path):
        self.bind_object.logger.info('start creating tables in sg_transcripts_step')
        transcript_id = ObjectId(transcript_id)
        step_files = os.listdir(statistics_path)
        data_list = list()
        for f in step_files:
            m = re.search(r'trans_count_stat_(\S+).txt', f)
            if m:
                step_list = list()
                step = m.group(1)
                files = os.path.join(statistics_path, f)
                fr = open(files)
                next(fr)
                for line in fr:
                    step_dic = dict()
                    step_range = line.strip().split('\t')[0]
                    num = line.strip().split('\t')[1]
                    step_dic[step_range] = num
                    step_list.append(step_dic)
                data = [
                    ('transcripts_id', transcript_id),
                    ('step', int(step)),
                    ('step_data', step_list)
                ]
                data = SON(data)
                data_list.append(data)
        self.create_db_table('sg_transcripts_step', data_list)
        self.bind_object.logger.info('succeed in creating tables in sg_transcripts_step')

    def add_transcripts_relations(self, transcript_id, statistics_path):
        self.bind_object.logger.info('start creating tables in sg_transcripts_relations')
        transcript_id = ObjectId(transcript_id)
        relation_files = os.listdir(statistics_path)
        dic = dict()
        data_list = list()
        name_list = list()
        for files in relation_files:
            m = re.search(r'(\S+)\.gtf\..+_([0-9]+)\.txt', files)
            if m:
                names = m.group(1)
                steps = m.group(2)
                if names not in name_list:
                    name_list.append(names)
                files_with_path = os.path.join(statistics_path, files)
                with open(files_with_path) as fr:
                    dic[(names, steps)] = list()
                    for line in fr:
                        gene_trans_dic = dict()
                        lines = line.strip().split('\t')
                        if len(lines) >= 2:
                            gene_trans_dic[lines[0]] = lines[1]
                            dic[(names, steps)].append(gene_trans_dic)
                    data = [
                        ('transcripts_id', transcript_id),
                        ('type_of_trans_or_genes', names),
                        ('step_of_trans_or_genes', int(steps)),
                        ('data_list', dic[(names, steps)]),
                    ]
                    data = SON(data)
                    data_list.append(data)
        self.create_db_table('sg_transcripts_relations', data_list)
        self.update_db_record('sg_transcripts', transcript_id, type_of_trans_or_genes=name_list)
        self.bind_object.logger.info('succeed in creating tables in sg_transcripts_relations')

    def add_transcripts_seq_type(self, transcript_id, code_file):
        self.bind_object.logger.info('start creating tables in sg_transcripts_seq_type')
        transcript_id = ObjectId(transcript_id)
        data_list = list()
        code_list = list()
        all_code = ['=', 'c', 'j', 'e', 'i', 'o', 'p', 'r', 'u', 'x', 's', '.']
        with open(code_file) as fr:
            for line in fr:
                new_gene_list = list()
                lines = line.strip().split('\t')
                gene_list = lines[1].strip().split(',')
                code_list.append(lines[0])
                for ids in gene_list:
                    if ids.startswith('transcript:'):
                        strinfo = re.compile('transcript:')
                        new_ids = strinfo.sub('', ids)
                        new_gene_list.append(new_ids)
                    else:
                        new_gene_list.append(ids)
                data = [
                    ('transcripts_id', transcript_id),
                    ('class_code', lines[0]),
                    ('num', int(lines[2])),
                    ('type', 'transcripts'),
                ]
                data = SON(data)
                data_list.append(data)
        for code in all_code:
            if code in code_list:
                pass
            else:
                data = [
                    ('transcripts_id', transcript_id),
                    ('class_code', code),
                    ('num', 0),
                    ('type', 'transcripts'),
                    ('gene_list', list()),
                ]
                data = SON(data)
                data_list.append(data)
        self.create_db_table('sg_transcripts_seq_type', data_list)
        self.update_db_record('sg_transcripts', transcript_id, seq_type=code_list)
        self.bind_object.logger.info('succeed in creating tables in sg_transcripts_seq_type')

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'assemble_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'lnc_rna.lnc_rna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = LncRnaTestApiWorkflow(wheet)
        wf.sheet.id = 'lnc_rna'
        wf.sheet.project_sn = 'lnc_rna'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('lnc_rna.assemble')
        params = {'method': 'stringtie', 'submit_location': 'transcripts', 'task_id': wf.sheet.id, 'task_type': 2}
        wf.test_api.add_assemble_result(
            params=params,
            all_gtf_path='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-stringtie/output/Stringtie',
            merged_path='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-stringtie/output/StringtieMerge',
            statistics_path='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-stringtie/output/Statistics'
        )

if __name__ == '__main__':
    unittest.main()
