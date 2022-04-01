# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from mainapp.controllers.project.small_rna_controller import SmallRnaController
import web
import re
import os
import json
from mainapp.libs.signature import check_sig
import sqlite3
import unittest

class QuerySeqAction(SmallRnaController):
    def __init__(self):
        super(QuerySeqAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()

        # check
        ret = self.check_params(data)
        if ret != True:
            return ret

        # download
        seq_db = self.get_seq_db(data)
        if not os.path.exists(seq_db):
            info = {'success': False, 'info': 'can not find local seq_db'}
            return json.dumps(info)

        # query
        cds, sequence, seq_len = self.query_seq(data.seq_type, data.seq_id, seq_db)

        # save
        mongo_data = dict([
            ('task_id', data.task_id),
            ('seq_id', data.seq_id),
            ('sequence', sequence),
            ('seq_len', seq_len),
            ('cds_info', cds),
            ('is_new', False)
        ])
        self.small_rna.insert_seq(mongo_data)
        info = {'success': True, 'info': 'query success'}
        return json.dumps(info)

    def check_params(self, data):
        expected_args = ['task_id', 'seq_id', 'seq_db', 'seq_type']
        for arg in expected_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)
        if data['seq_type'].upper() not in ['T', 'G']:
            info = {'success': False, 'info': 'The value of seq_type must be T or G'}
            return json.dumps(info)
        m = re.match(r'^\w+://\S+/.+$', data['seq_db'])
        if not m:
            info = {'success': False, 'info': 'The value of seq_db must be s3 remote input path'}
            return json.dumps(info)
        return True

    def get_seq_db(self, data):
        inter_dir = self.create_tmp_dir(data['task_id'], 'query_seq/')
        download_file = self.download_from_s3(data['seq_db'], inter_dir)
        return download_file

    @staticmethod
    def query_seq(seq_type, seq_id, seq_db):
        '''
        if seq_type is G, cds_data is None, return only the sequence and length of sequence
        if seq_type is T, cds_data contains the cds and pep information of related transcript
        '''
        conn = sqlite3.connect(seq_db)
        cursor = conn.cursor()

        # get cds and pep sequence
        if seq_type.upper() == 'T':
            cursor.execute("SELECT * FROM {} WHERE {}='{}'".format('trans_annot', 'transcript_id', seq_id))
            # type of cds_info is list
            # [(transcript_id, cds_id, pep_id, cds_seq, pep_seq, orf_type), (...), ...]
            cds_info = cursor.fetchall()
            cds_data = list()
            fileds = ['transcript_id', 'cds_id', 'pep_id', 'cds_seq', 'pep_seq', 'orf_type']
            if cds_info:
                # each is a tuple which contains six items correspond to the six elements in fileds
                for each in cds_info:
                    tmp = dict(zip(fileds, each))
                    if tmp['pep_seq'] == 'None':
                        tmp['pep_length'] = 0
                    else:
                        tmp['pep_length'] = len(tmp['pep_seq'])
                    tmp['cds_length'] = len(tmp['cds_seq'])
                    cds_data.append(tmp)
            else:
                print '{} not found'.format(seq_id)
                tmp = dict(zip(fileds, [seq_id, 'None', 'None', 'None', 'None', 'None',]))
                tmp['cds_length'] = 0
                tmp['pep_length'] = 0
                cds_data.append(tmp)
        else:
            cds_data = None

        # get gene or transcript sequence
        if seq_type.upper() == 'G':
            query_table = 'gene_seq'
        elif seq_type.upper() == 'T':
            query_table = 'trans_seq'
        try:
            cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'seq_id', seq_id))
        except:
            # 无参数据库表不一样
            query_table = 'transcript_gene'
            cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'g_id', seq_id))
        seq = cursor.fetchall()
        if seq:
            sequence = seq[0][0]
        else:
            sequence = 'None'
        cursor.close()

        return cds_data, sequence, len(sequence)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test_G(self):
        import os
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 'i/small_rna/query_seq '
        cmd += '-b http://192.168.12.102:9090 '
        args = dict(
            task_id='tsg_33087',
            seq_id='ENSG00000282222',
            seq_db='s3://smallrna/files/m_188/188_5c08b2758d630/tsg_33087/workflow_results/Sequence_database/refrna_seqs.db',
            seq_type='G',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

    def test_T(self):
        import os
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 'i/small_rna/query_seq '
        cmd += '-b http://192.168.12.102:9090 '
        args = dict(
            task_id='tsg_33087',
            seq_id='ENST00000349238',
            seq_db='s3://smallrna/files/m_188/188_5c08b2758d630/tsg_33087/workflow_results/Sequence_database/refrna_seqs.db',
            seq_type='T',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
