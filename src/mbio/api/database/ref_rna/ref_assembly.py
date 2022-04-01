# -*- coding: utf-8 -*-
# __author__ = 'zhaoyue.wang'
import json
import datetime
import os
import re
from biocluster.api.database.base import Base, report_check
from bson.son import SON
import types
from bson.objectid import ObjectId
import pymongo
from biocluster.config import Config


class RefAssembly(Base):
    def __init__(self, bind_object):
        super(RefAssembly, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_assembly_result(self, name=None, params=None, all_gtf_path=None, merged_path=None, Statistics_path=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        merged_list = []
        merged_gtf_path = dict()
        merged_gtf_path["gtf"] = merged_path + '/merged.gtf'
        merged_list.append(merged_gtf_path)
        merged_fa_path = dict()
        merged_fa_path["fa"] = merged_path + '/merged.fa'
        merged_list.append(merged_fa_path)
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'Assembly_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'end',
            'desc': 'cufflinks拼接组装结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'step': [200, 300, 600, 1000],
            'seq_type': ["i", "j", "o", "u", "x"],
            'merge_path': merged_list,
            'Sample_gtf_path': all_gtf_path,
            "type_of_trans_or_genes": ["old_gene", "old_trans", "new_gene", "new_trans"],
            "step_of_trans_or_genes": [1, 5, 10, 20],

        }
        collection = self.db['sg_transcripts']
        transcript_id = collection.insert_one(insert_data).inserted_id
        code_files = Statistics_path + '/code_num.txt'
        self.add_transcripts_step(transcript_id=transcript_id, Statistics_path=Statistics_path)
        self.add_transcripts_relations(transcript_id=transcript_id, Statistics_path=Statistics_path)
        self.add_transcripts_seq_type(transcript_id=transcript_id, code_file=code_files)
        return transcript_id

    def add_transcripts_step(self, transcript_id, Statistics_path):

        if not isinstance(transcript_id, ObjectId):
            if isinstance(transcript_id, types.StringTypes):
                transcript_id = ObjectId(transcript_id)
            else:
                raise Exception('transcript_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(Statistics_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(Statistics_path))

        step_files = os.listdir(Statistics_path)
        data_list = []
        for f in step_files:
            m = re.search(r'trans_count_stat_(\S+)\.txt', f)
            if m:
                step_list = []
                step = m.group(1)
                files = os.path.join(Statistics_path, f)
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
        try:
            collection = self.db['sg_transcripts_step']
            collection.insert_many(data_list)
        except:
            raise Exception("导入步长信息失败!")
        else:
            self.bind_object.logger.info("导入步长信息成功!")

    def add_transcripts_relations(self, transcript_id, Statistics_path):

        if not isinstance(transcript_id, ObjectId):
            if isinstance(transcript_id, types.StringTypes):
                transcript_id = ObjectId(transcript_id)
            else:
                raise Exception('transcript_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(Statistics_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(Statistics_path))
        relation_files = os.listdir(Statistics_path)
        dic = {}
        data_list = []
        name_list = []
        for files in relation_files:
            m = re.search(r'(\S+)\.gtf\..+_([0-9]+)\.txt', files)
            if m:
                names = m.group(1)
                steps = m.group(2)
                name_list.append(names)
                files_with_path = Statistics_path + "/" + files
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
        try:
            collection = self.db['sg_transcripts_relations']
            collection.insert_many(data_list)
            collection = self.db['sg_transcripts']
            collection.update({'_id': ObjectId(transcript_id)}, {'$set': {'type_of_trans_or_genes': name_list}})
        except:
            raise Exception("导入class_code信息：失败!")
        else:
            self.bind_object.logger.info("导入class_code信息：成功!")

    def add_transcripts_seq_type(self, transcript_id, code_file):
        if not isinstance(transcript_id, ObjectId):
            if isinstance(transcript_id, types.StringTypes):
                transcript_id = ObjectId(transcript_id)
            else:
                raise Exception('transcript_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(code_file):
            raise Exception('{}所指定的文件不存在，请检查！'.format(code_file))
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
                    ('gene_list', new_gene_list),
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
        try:
            collection = self.db['sg_transcripts_seq_type']
            collection.insert_many(data_list)
            collection = self.db['sg_transcripts']
            collection.update({'_id': ObjectId(transcript_id)}, {'$set': {'seq_type': code_list}})
        except:
            raise Exception("导入class_code信息：%s失败!")
        else:
            self.bind_object.logger.info("导入class_code信息成功!")

if __name__ == "__main__":

    a = RefAssembly()
    all_gtf_path = '/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_cufflinks_true_file/RefrnaAssemble/output/Cufflinks'
    merged_path = "/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_cufflinks_true_file/RefrnaAssemble/output/Cuffmerge"
    transcript_id = a.add_assembly_result(all_gtf_path=all_gtf_path, merged_path=merged_path)
    Statistics_path = '/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_cufflinks_true_file/RefrnaAssemble/output/Statistics'
    a.add_transcripts_step(transcript_id=transcript_id, Statistics_path=Statistics_path)
    a.add_transcripts_relations(transcript_id=transcript_id, Statistics_path=Statistics_path)
    code_files = '/mnt/ilustre/users/sanger-dev/workspace/20170508/Single_assembly_module_tophat_cufflinks_true_file/RefrnaAssemble/output/Statistics/code_num.txt'
    a.add_transcripts_seq_type(transcript_id=transcript_id, code_file=code_files)
