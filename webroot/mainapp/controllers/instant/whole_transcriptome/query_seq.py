# -*- coding: utf-8 -*-
import os
import web
import json
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig
import sqlite3
# __author__ = 'gdq'
import unittest
import re
import datetime

class QuerySeqAction(WholeTranscriptomeController):
    def __init__(self):
        super(QuerySeqAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        # check
        for arg in ['seq_id', 'category', 'level', 'task_id']:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': 'Argument:%s required' % arg, 'code': 'C2902101', 'variables': variables}
                return json.dumps(info)

        if data['level'] not in ['T', 'G']:
            info = {'success': False, 'info': 'unexpected level', 'code': 'C2902102', 'variables': ''}
            return json.dumps(info)

        genes_dict = self.whole_transcriptome.get_main_info_by_record('genes', task_id=data.task_id)
        # print "gene_dict is {}".format(genes_dict)
        task_info = self.whole_transcriptome.get_task_info(data.task_id)

        if os.path.exists(genes_dict['refrna_seqdb']):
            seq_dir = genes_dict['refrna_seqdb']
        else:
            seq_dir = task_info['output'] + '/other/gene_detail'

        if data['category'] == "mRNA":
            self.db_dir  = seq_dir + '/mrna/mrna.db'
        elif data['category'] == "lncRNA":
            self.db_dir  = seq_dir + '/lncrna/lncrna.db'
        elif data['category'] == "circRNA":
            self.db_dir  = seq_dir + '/lncrna/circrna.db'

        if re.match(r'^\w+://\S+/.+$', self.db_dir):
            inter_dir = self.create_tmp_dir(data.task_id, "query_seq/")
            download_seq = self.download_from_s3(self.db_dir, inter_dir=inter_dir)
            self.db_dir = os.path.join(download_seq)

        elif not os.path.exists(self.db_dir):
            variables = []
            variables.append(self.db_dir)
            info = {'success': False, 'info': '%s +  not exist' % self.db_dir, 'code': 'C2902103', 'variables': variables }
            return json.dumps(info)
        # print "data is {}".format(datetime.datetime.now())

        # query
        cds, sequence, seq_len = self.query_seq(data.level, data.seq_id, self.db_dir)
        # print "data is {}".format(datetime.datetime.now())
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
        self.whole_transcriptome.insert_seq(mongo_data)
        print "data is {}".format(datetime.datetime.now())
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
            cursor.execute("SELECT * FROM {} WHERE {}='{}'".format('sub_info', 'transcript_id', seq_id))
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
        query_table = 'gene_seq' if level.upper() == 'G' else 'transcript_seq'
        cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'seq_id', seq_id))
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
        cmd += "i/whole_transcriptome/query_seq "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="whole_transcriptome",
            seq_id="ENST00000642491",
            category="mRNA",
            level="T",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
