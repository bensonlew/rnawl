# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue,shicaiping'

import datetime
import os
import re
import types
import unittest

from bson.objectid import ObjectId
from bson.son import SON

from biocluster.api.database.base import report_check
from biocluster.config import Config
from mbio.api.database.ref_rna_v2.api_base import ApiBase


class RefAssembly(ApiBase):
    def __init__(self, bind_object):
        super(RefAssembly, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'
        _ = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]

    @report_check
    def add_assembly_result(self, params, all_gtf_path, merged_path, statistics_path):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        self.bind_object.logger.info('gtf store place: {}'.format(all_gtf_path))
        self.bind_object.logger.info('task id: {}'.format(task_id))
        self.bind_object.logger.info('project sn: {}'.format(project_sn))
        merged_list = []
        merged_gtf_path = dict()
        merged_gtf_path["gtf"] = merged_path + '/merged.gtf'
        merged_list.append(merged_gtf_path)
        merged_fa_path = dict()
        merged_fa_path["fa"] = merged_path + '/merged.fa'
        merged_list.append(merged_fa_path)
        name = 'Assembly_{}'.format(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'params': params,
            'status': 'end',
            'desc': 'cufflinks拼接组装结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'step': [200, 300, 600, 1000],
            'seq_type': ["i", "j", "o", "u", "x"],
            "step_of_trans_or_genes": [1, 5, 10, 20],

        }
        collection = self.db['sg_transcripts']
        transcript_id = collection.insert_one(insert_data).inserted_id
        code_files = statistics_path + '/code_num.txt'
        self.add_transcripts_step(transcript_id=transcript_id, statistics_path=statistics_path)
        self.add_transcripts_relations(transcript_id=transcript_id, statistics_path=statistics_path)
        self.add_transcripts_seq_type(transcript_id=transcript_id, code_file=code_files)
        self.update_db_record('sg_transcripts', transcript_id, status="end", main_id=transcript_id)
        return transcript_id

    def add_transcripts_step(self, transcript_id, statistics_path):
        if not isinstance(transcript_id, ObjectId):
            if isinstance(transcript_id, types.StringTypes):
                transcript_id = ObjectId(transcript_id)
            else:
                self.bind_object.set_error('transcript_id必须为ObjectId对象或其对应的字符串！', code="53701910")
        if not os.path.exists(statistics_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=statistics_path, code="53701911")
        step_files = os.listdir(statistics_path)
        data_list = []
        for f in step_files:
            m = re.search(r'trans_count_stat_(\S+)\.txt', f)
            if m:
                step_list = []
                step = m.group(1)
                files = os.path.join(statistics_path, f)
                fr = open(files, "r")
                next(fr)
                for line in fr:
                    step_dic = dict()
                    step_range = line.strip().split("\t")[0]
                    num = line.strip().split("\t")[1]
                    step_dic[step_range] = num
                    step_list.append(step_dic)
                data = [
                    ('transcripts_id', transcript_id),
                    ('step', int(step)),
                    ('step_data', step_list)
                ]
                data = SON(data)
                data_list.append(data)
        collection = self.db['sg_transcripts_step']
        collection.insert_many(data_list)

    def add_transcripts_relations(self, transcript_id, statistics_path):
        if not isinstance(transcript_id, ObjectId):
            if isinstance(transcript_id, types.StringTypes):
                transcript_id = ObjectId(transcript_id)
            else:
                self.bind_object.set_error('transcript_id必须为ObjectId对象或其对应的字符串！', code="53701913")
        if not os.path.exists(statistics_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(statistics_path), code="53701914")
        relation_files = os.listdir(statistics_path)
        dic = {}
        data_list = []
        name_list = []
        for files in relation_files:
            m = re.search(r'(\S+)\.gtf\..+_([0-9]+)\.txt', files)
            if m:
                names = m.group(1)
                steps = m.group(2)
                if names not in name_list:
                    name_list.append(names)
                files_with_path = statistics_path + "/" + files
                with open(files_with_path, "r") as fr:
                    dic[(names, steps)] = []
                    for line in fr:
                        gene_trans_dic = dict()
                        lines = line.strip().split("\t")
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
        collection = self.db['sg_transcripts_relations']
        collection.insert_many(data_list)
        collection = self.db['sg_transcripts']
        collection.update({'_id': ObjectId(transcript_id)}, {'$set': {'type_of_trans_or_genes': name_list}})

    def add_transcripts_seq_type(self, transcript_id, code_file):
        if not isinstance(transcript_id, ObjectId):
            if isinstance(transcript_id, types.StringTypes):
                transcript_id = ObjectId(transcript_id)
            else:
                self.bind_object.set_error('transcript_id必须为ObjectId对象或其对应的字符串！', code="53701916")
        if not os.path.exists(code_file):
            self.bind_object.set_error('%s所指定的文件不存在，请检查！', variables=code_file, code="53701917")
        data_list = []
        code_list = []
        all_code = ['=', 'c', 'j', 'e', 'i', 'o', 'p', 'r', 'u', 'x', 's', '.']
        with open(code_file, "r") as fr:
            for line in fr:
                new_gene_list = []
                lines = line.strip().split("\t")
                gene_list = lines[1].strip().split(",")
                code_list.append(lines[0])
                for ids in gene_list:
                    if ids.startswith("transcript:"):
                        strinfo = re.compile('transcript:')
                        new_ids = strinfo.sub('', ids)
                        new_gene_list.append(new_ids)
                    else:
                        new_gene_list.append(ids)
                data = [
                    ('transcripts_id', transcript_id),
                    ('class_code', lines[0]),
                    ('num', int(lines[2])),
                    ('type', "transcripts"),
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
                    ('type', "transcripts"),
                    ('gene_list', []),
                ]
                data = SON(data)
                data_list.append(data)
        collection = self.db['sg_transcripts_seq_type']
        collection.insert_many(data_list)
        collection = self.db['sg_transcripts']
        collection.update({'_id': ObjectId(transcript_id)}, {'$set': {'seq_type': code_list}})


class TestFunction(unittest.TestCase):
    """
    This is test for the api. Just run this script to do test.
    """

    def test(self):
        import json
        import random
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet

        task_id = 'i-sanger_218717'
        project_sn = '35910_5dc8f4043c529'

        db = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
        main_doc = db['sg_transcripts'].find_one({'task_id': task_id})
        if main_doc and 'main_id' in main_doc:
            transcripts_id = main_doc['main_id']
            db['sg_transcripts_relations'].delete_many({'transcripts_id': transcripts_id})
            db['sg_transcripts_seq_type'].delete_many({'transcripts_id': transcripts_id})
            db['sg_transcripts_step'].delete_many({'transcripts_id': transcripts_id})
            db['sg_transcripts'].delete_many({'task_id': task_id})

        data = {
            'id': 'ref_assembly_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v2.refrna_test_api',
            'options': {}
        }
        params = json.dumps({'task_id': task_id, 'submit_location': 'transcripts', 'task_type': 2}, sort_keys=True)
        refrna_assemble_output_dir = '/mnt/lustre/users/sanger/workspace/20200304/Refrna_majorbio_242245' \
                                     '/RefrnaAssemble/output'
        all_gtf_path = os.path.join(refrna_assemble_output_dir, 'Stringtie')
        merged_path = os.path.join(refrna_assemble_output_dir, 'StringtieMerge')
        statistics_path = os.path.join(refrna_assemble_output_dir, 'Statistics')
        wheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.sheet.id = task_id
        wf.sheet.project_sn = project_sn
        wf.test_api = wf.api.api('ref_rna_v2.ref_assembly')
        wf.test_api.add_assembly_result(params=params, all_gtf_path=all_gtf_path, merged_path=merged_path,
                                        statistics_path=statistics_path)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
