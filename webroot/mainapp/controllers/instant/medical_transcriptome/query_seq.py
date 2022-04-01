# -*- coding: utf-8 -*-
import os
import web
import json
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
import sqlite3
# __author__ = 'gdq'
import unittest
import re

class QuerySeqAction(MedicalTranscriptomeController):
    def __init__(self):
        super(QuerySeqAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        # check
        for arg in ['seq_id', 'seq_db', 'level', 'task_id']:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': 'Argument:%s required' % arg, 'code': 'C2902101', 'variables': variables}
                return json.dumps(info)

        if data['level'] not in ['T', 'G']:
            info = {'success': False, 'info': 'unexpected level', 'code': 'C2902102', 'variables': ''}
            return json.dumps(info)

        if re.match(r'^\w+://\S+/.+$', data['seq_db']) or re.match(r'/mnt/ilustre', data['seq_db']):
            inter_dir = self.create_tmp_dir(data.task_id, "query_seq/")
            download_seq = self.download_from_s3(data['seq_db'], inter_dir=inter_dir)
            data['seq_db'] = os.path.join(download_seq)

        elif not os.path.exists(data['seq_db']):
            variables = []
            variables.append(data['seq_db'])
            info = {'success': False, 'info': '%s +  not exist' % data['seq_db'], 'code': 'C2902103', 'variables': variables }
            return json.dumps(info)

        print("lalala")

        # query
        cds, sequence, seq_len = self.query_seq(data.level, data.seq_id, data.seq_db)
        if data.seq_id.startswith('MSTR') or data.seq_id.startswith('TCON') or data.seq_id.startswith('XLOC') :
            is_new = True
        else:
            is_new = False
        # save
        mongo_data = dict([('task_id', data.task_id),
                           ('seq_id', data.seq_id),
                           ('sequence', sequence),
                           ('seq_len', seq_len),
                           ('cds_info', cds),
                           ('is_new', is_new)
                           ])
        self.medical_transcriptome.insert_seq(mongo_data)
        info = {"success": True, "info": "query success", 'code': 'C2902104', 'variables': ''}
        return json.dumps(info)

    @staticmethod
    def query_seq(level, seq_id, seq_db):
        """
        当查询的是基因时，则仅仅返回基因的序列信息
        当查询的是转录本时，则返回转录本序列信息，并且返回该转录本对应的所有cds和pep信息
        """
        conn = sqlite3.connect(seq_db)
        cursor = conn.cursor()

        # get cds and pep sequence
        if level.upper() == 'T':
            cursor.execute("SELECT * FROM {} WHERE {}='{}'".format('trans_annot', 'transcript_id', seq_id))
            cds_info = cursor.fetchall()
            cds_data = list()
            fileds = ['transcript_id', 'cds_id', 'pep_id', 'cds_seq', 'pep_seq', 'orf_type']
            if cds_info:
                for each in cds_info:
                    tmp = dict(zip(fileds, each))
                    if tmp['pep_seq'] == 'None':
                        tmp['pep_length'] = 0
                    else:
                        tmp['pep_length'] = len(tmp['pep_seq'])
                    tmp['cds_length'] = len(tmp['cds_seq'])
                    cds_data.append(tmp)
            else:
                print("{} not found".format(seq_id))
                tmp = dict(zip(fileds, [seq_id, 'None', 'None', 'None', 'None', 'None',]))
                tmp['cds_length'] = 0
                tmp['pep_length'] = 0
                cds_data.append(tmp)
        else:
            cds_data = None

        # get gene or transcript sequence
        query_table = 'gene_seq' if level.upper() == 'G' else 'trans_seq'
        cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'seq_id', seq_id))
        seq = cursor.fetchall()
        if seq:
            sequence = seq[0][0]
        else:
            ### 部分项目之前分隔符写错造成序列id出错
            cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'seq_id', seq_id + "("))
            seq = cursor.fetchall()
            if seq:
                sequence = seq[0][0]
            else:
                sequence = 'None'
        cursor.close()

        return cds_data, sequence, len(sequence)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        import os
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/medical_transcriptome/query_seq "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="gbnijnfcuslii6rieegpbdgpt5",
            seq_id="ENSG00000282222",
            seq_db="s3://common/files/m_188/jf9c6bo3ql248thi1taa4c4plu/gbnijnfcuslii6rieegpbdgpt5/intermediate_results/SequenceDatabase/refrna_seqs.db",
            level="G",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
