# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.protein_transcript_dia_v3.api_base import ApiBase
import json
from bson.objectid import ObjectId
import types
import sys
from biocluster.config import Config


# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class ProteinTranscript(ApiBase):
    def __init__(self, bind_object):
        super(ProteinTranscript, self).__init__(bind_object)
        self._project_type = 'dia'

    #@report_check
    def add_relation(self, main_id = None, params=None, project_sn='protein_transcript_dia_v3',
                    task_id='protein_transcript_dia_v3', relation_file=None, table = 'sg_p2g_relationship', detail_id = 'rela_id'):
        insert_data = dict()
        rela_data = dict()
        with open(relation_file, 'r') as r_f:
            for line in r_f.readlines():
                line = line.strip().split('\t')
                if line:
                    try:
                        rela_data[line[0]] = line[1]
                        rela_data['{}_num'.format(line[0])] = len(line[1].split(';'))
                    except:
                        rela_data[line[0]] = ''
                        rela_data['{}_num'.format(line[0])] = 0
        if params is None:
            params = {"software": "blastp"}
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': self.bind_object.sheet.id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'protein_transcript_dia_v3_relationship',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            }
        if not main_id:
            # insert_data.update(rela_data)
            rel_id = self.create_db_table(table, [insert_data])
            rela_data.update({detail_id : rel_id})
            self.update_db_record(table+'_prvenn', rel_id, status="end", main_id=rel_id)
        else:
            if not isinstance(main_id, ObjectId):
                if isinstance(main_id, types.StringTypes):
                    main_id = ObjectId(main_id)
                else:
                    raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
            # self.update_db_record(table, main_id, **rela_data)
            self.update_db_record(table, main_id, status="end", main_id=main_id)
            rela_data.update({detail_id: main_id})
            # if sys.getsizeof(rela_data) < 16793598:
            #     self.create_db_table(table+'_prvenn', [rela_data])
            # else:
            #     for k, v in rela_data.items():
            #         if k == 'failed_proteins':
            #             rel_id = self.create_db_table(table + '_prvenn', [{k,v}])
            #     for k, v in rela_data.items():
            #         if k != 'failed_proteins':
            #             tmp = {k:v}
            #             print(tmp)
            #             self.update_db_record(table + '_prvenn', rel_id, **tmp)

            # fix venn gene num showed is different from transcriptome data ,which caused by one gene correspond to multiple proteins. By xuxi on 20211008.
            try:
                gene_list_tmp = list(map(list, zip(*[i.split('|') for i in rela_data["related"].split(';')])))[1]
                rela_data["failed_transcripts_num"] -= (len(gene_list_tmp)-len(set(gene_list_tmp)))
            except:
                pass
            #

            for k, v in rela_data.items():
                if k == detail_id:
                    tmp = {k: v}
                    rel_id = self.create_db_table(table + '_prvenn', [tmp])
                    print(rel_id)
            for k, v in rela_data.items():
                if k != detail_id:
                    try:
                        num = len(v.split(';'))
                    except:
                        num = 0
                    if num > 15000:
                        v = 'too_large'
                    tmp = {k: v}
                    print(tmp)
                    self.update_db_record(table + '_prvenn', rel_id, **tmp)

    def update_task_relation(self, main_id):
        main_info = self.get_main_record_by_kwargs('sg_p2g_relationship', main_id=ObjectId(main_id))
        data = ({
            'task_id': main_info['task_id'],
            'project_sn': main_info['project_sn'],
            'project_type': 'dia',
            'relation_task_id': main_info['rna_task_id'],
            'relation_project_type': main_info['rna_type'],
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'created_id': 'proteni_transcript_dia_v3_controller',
            'is_delete': 0})
        db = Config().get_mongo_client(mtype='project')[Config().get_mongo_dbname(mtype='project')]
        my_collection = db['sg_task_relations']
        my_collection.insert_one(data)



