# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.protein_transcript.api_base import ApiBase
import json
from bson.objectid import ObjectId
import types
import sys


# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class ProteinTranscript(ApiBase):
    def __init__(self, bind_object):
        super(ProteinTranscript, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    #@report_check
    def add_relation(self, main_id = None, params=None, project_sn='protein_transcript',
                    task_id='protein_transcript', relation_file=None, table = 'sg_p2g_relationship', detail_id = 'rela_id'):

        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
        rela_dict = dict()
        with open(relation_file, 'r') as r_f:
            for line in r_f.readlines():
                if line:
                    try:
                        category, name = line.strip().split('\t')
                        rela_dict[category] = name
                        rela_dict['{}_num'.format(category)] = len(name.split(','))
                    except:
                        category = line.strip().split('\t')[0]
                        rela_dict[category] = ''
                        rela_dict['{}_num'.format(category)] = 0

        gene_str = rela_dict['failed_transcripts'] + ',' + rela_dict['related']
        gene_str = gene_str.strip(',')
        protein_str = rela_dict['failed_proteins'] + ',' + rela_dict['related']
        protein_str = protein_str.strip(',')
        data_list = list()
        data_list.append({'names': gene_str, 'category': 'gene', 'type': 'venn', detail_id: main_id})
        data_list.append({'names': protein_str, 'category': 'protein', 'type': 'venn', detail_id: main_id})
        self.create_db_table(table + '_prvenn', data_list)
        venn_data = {'names': 'names', 'category': 'category', 'condition': {'type': 'venn'}}
        venn_data_json = json.dumps(venn_data, sort_keys=True, separators=(',', ':'))
        rela_list = ['failed_transcripts_num', 'failed_proteins_num', 'related_num']
        data_dict2 = {detail_id: main_id}
        for i in rela_list:
            data_dict2.update({i: rela_dict[i]})
        data_list2 = [data_dict2]
        print data_list2
        columns_list = list()
        columns_list.append({'field': 'failed_transcripts_num', 'filter': False, 'sort': False, 'title': 'Gene number', 'type': 'int'})
        columns_list.append({'field': 'failed_proteins_num', 'filter': False, 'sort': False, 'title': 'Protein number', 'type': 'int'})
        columns_list.append({'field': 'related_num', 'filter': False, 'sort': False, 'title': 'Related number', 'type': 'int'})
        data_columns = {'column': columns_list, 'condition': {}}
        columns_data = json.dumps(data_columns)
        self.create_db_table(table + '_detail', data_list2)
        self.update_db_record(table, main_id, status="end", main_id=main_id, venn_data=venn_data_json, column_data_detail=columns_data)




        #
        #
        #
        # rela_data = dict()
        # with open(relation_file, 'r') as r_f:
        #     for line in r_f.readlines():
        #         line = line.strip().split('\t')
        #         if line:
        #             try:
        #                 rela_data[line[0]] = line[1]
        #                 rela_data['{}_num'.format(line[0])] = len(line[1].split(';'))
        #             except:
        #                 rela_data[line[0]] = ''
        #                 rela_data['{}_num'.format(line[0])] = 0
        # if params is None:
        #     params = {"software": "blastp"}
        # params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        # insert_data = {
        #     'project_sn': self.bind_object.sheet.project_sn,
        #     'task_id': self.bind_object.sheet.id,
        #     'params': params if params else "",
        #     'status': 'start',
        #     'desc': 'protein_transcript_relationship',
        #     'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
        #     }
        # if not main_id:
        #
        #     main_id = self.create_db_table(table, [insert_data])
        #     print 'main_ID'
        #     print main_id
        #
        # if not isinstance(main_id, ObjectId):
        #     if isinstance(main_id, types.StringTypes):
        #         main_id = ObjectId(main_id)
        #     else:
        #         raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
        #
        # rela_data.update({detail_id: main_id})
        # print 'rela_data'
        # print rela_data
        #
        # for k, v in rela_data.items():
        #     if k == 'failed_proteins':
        #         tmp = {k: v}
        #         rel_id = self.create_db_table(table + '_prvenn', [tmp])
        #         print 'rel_ID'
        #         print(rel_id)
        # for k, v in rela_data.items():
        #     if k != 'failed_proteins':
        #         try:
        #             num = len(v.split(';'))
        #         except:
        #             num = 0
        #         if num > 15000:
        #             v = 'too_large'
        #         tmp = {k: v}
        #         print 'TMP'
        #         print(tmp)
        #         self.update_db_record(table + '_prvenn', rel_id, **tmp)
        # self.update_db_record(table, main_id, status="end", main_id=main_id)



