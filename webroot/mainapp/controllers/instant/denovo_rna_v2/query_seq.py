# -*- coding: utf-8 -*-
# __author__ = 'gdq'

import os
import web
import json
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mainapp.libs.signature import check_sig
import sqlite3
import re

class QuerySeqAction(DenovoRnaV2Controller):
    def __init__(self):
        super(QuerySeqAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        # check
        for arg in ['seq_id', 'seq_db', 'seq_type', 'task_id']:
            if not hasattr(data, arg):
                var = []
                var.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % (arg), "code": 'C1601801', "variables": var}
                return json.dumps(info)

        if data['seq_type'] not in ['T', 'G']:
            info = {'success': False, 'info': 'unexpected seq_type', "code": 'C1601802', "variables": ''}
            return json.dumps(info)

        if re.match(r'^\w+://\S+/.+$', data['seq_db']) or re.match(r'/mnt/ilustre', data['seq_db']):
            inter_dir = self.create_tmp_dir(data.task_id, 'query_seq/')
            download_seq = self.download_from_s3(data['seq_db'], inter_dir=inter_dir)
            data['seq_db'] = os.path.join(download_seq)
        elif not os.path.exists(data['seq_db']):
            var = []
            var.append(data['seq_db'])
            info = {'success': False, 'info': data['seq_db'] + ' not exist', "code": 'C1601803', "variables": var}
            return json.dumps(info)

        # query
        cds, sequence, seq_len = self.query_seq(data.seq_type, data.seq_id, data.seq_db)
        # save
        mongo_data = dict([('task_id', data.task_id),
                           ('seq_id', data.seq_id),
                           ('sequence', sequence),
                           ('seq_len', seq_len),
                           ('cds_info', cds),
                           ])
        self.denovo_rna_v2.insert_seq(mongo_data)
        info = {"success": True, "info": "query success"}
        return json.dumps(info)

    @staticmethod
    def query_seq(seq_type, seq_id, seq_db):
        conn = sqlite3.connect(seq_db)
        cursor = conn.cursor()
        if seq_type == 'T':
            query_field = 't_id'
        else:
            query_field = 'g_id'

        cursor.execute("SELECT * FROM {} WHERE {}='{}'".format('seq_cds', query_field, seq_id))
        cds_info = cursor.fetchall()
        cds_data = list()
        if cds_info:
            for each in cds_info:
                tmp = dict(zip(['t_id', 'g_id', 'cds_pos', 'cds', 'pep', 'type'], each))
                tmp['cds_length'] = len(tmp['cds'])
                tmp['pep_length'] = len(tmp['pep'])
                cds_data.append(tmp)
        else:
            print("{} not found".format(seq_id))
            tmp = dict(zip(['t_id', 'g_id', 'cds_pos', 'cds', 'pep', 'type'],
                           [seq_id, seq_id, 'None', 'None', 'None', 'None']))
            tmp['cds_length'] = 0
            tmp['pep_length'] = 0
            cds_data.append(tmp)
        #
        cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format('transcript_gene', query_field, seq_id))
        seq = cursor.fetchall()
        if seq:
            sequence = seq[0][0]
        else:
            sequence = 'None'
        cursor.close()
        return cds_data, sequence, len(sequence)
