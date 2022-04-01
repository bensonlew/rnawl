# -*- coding: utf-8 -*-
import os
import web
import json
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mainapp.libs.signature import check_sig
import sqlite3
import unittest
# __author__ = 'guhaidong'


class MetabsetAction(ItraqTmtController):
    def __init__(self):
        super(MetabsetAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        for arg in ['seq_id', 'seq_db', 'task_id']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)

        if not os.path.exists(data['seq_db']):
            info = {'success': False, 'info': data['seq_db'] + ' not exist'}
            return json.dumps(info)

        # query
        sequence, seq_len = self.query_seq(data.seq_id, data.seq_db)
        # save
        mongo_data = dict([('task_id', data.task_id),
                           ('seq_id', data.seq_id),
                           ('sequence', sequence),
                           ('seq_len', seq_len),
                           ])
        self.itraq_tmt.insert_seq(mongo_data)
        info = {"success": True, "info": "query success", }
        return json.dumps(info)

    @staticmethod
    def query_seq(seq_id, seq_db):
        conn = sqlite3.connect(seq_db)
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM {} WHERE pep_id ='{}'".format('seq_protein', seq_id))
        seq = cursor.fetchall()
        if seq:
            sequence = seq[0][1]
        else:
            sequence = 'None'
        cursor.close()
        return sequence, len(sequence)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/itraq_and_tmt/query_seq "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="itraq1",
            seq_id="K4BQC4",
            seq_db="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/itraq_test/seq_db.sqlite3",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
