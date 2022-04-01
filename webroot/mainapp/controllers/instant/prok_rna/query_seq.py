# -*- coding: utf-8 -*-
import os
import web
import json
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mainapp.libs.signature import check_sig
import sqlite3
import unittest
import re
# __author__ = 'scp'


class QuerySeqAction(ProkRNAController):
    def __init__(self):
        super(QuerySeqAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        # check
        for arg in ['seq_id', 'task_id']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)

        rock_index = self.prok_rna.get_task_info(data.task_id)['rock_index']
        if re.match(r'^\w+://\S+/.+$', rock_index):
            inter_dir = self.create_tmp_dir(data.task_id, "rock_index/")
            rock_index_fa = self.download_from_s3(os.path.join(rock_index, 'cds.fa.db.sqlite3'), inter_dir=inter_dir)
            rock_index_faa = self.download_from_s3(os.path.join(rock_index, 'cds.faa.db.sqlite3'), inter_dir=inter_dir)
            rock_index = rock_index_faa.split('cds.faa.db.sqlite3')[0]

        if not os.path.exists(os.path.join(rock_index, 'cds.fa.db.sqlite3')):
            info = {'success': False, 'info': '基因db不存在'}
            return json.dumps(info)
        if not os.path.exists(os.path.join(rock_index, 'cds.faa.db.sqlite3')):
            info = {'success': False, 'info': '蛋白db不存在'}
            return json.dumps(info)

        # query
        gene_seq, gene_len, pro_seq, pro_len = self.query_seq(data.seq_id, rock_index)
        # save
        mongo_data = dict([('task_id', data.task_id),
                           ('seq_id', data.seq_id),
                           ('gene_seq', gene_seq),
                           ('gene_len', gene_len),
                           ('protein_seq', pro_seq),
                           ('protein_len', str(int(pro_len) - 1)),
                           ])
        self.prok_rna.insert_seq(mongo_data)
        info = {"success": True, "info": "query success"}
        return json.dumps(info)

    @staticmethod
    def query_seq(seq_id, db_path):
        conn_gene = sqlite3.connect(os.path.join(db_path, 'cds.fa.db.sqlite3'))
        conn_protein = sqlite3.connect(os.path.join(db_path, 'cds.faa.db.sqlite3'))
        cursor_gene = conn_gene.cursor()
        cursor_protein = conn_protein.cursor()
        cursor_gene.execute("SELECT * FROM {} WHERE gene_id ='{}'".format('seq_gene', seq_id))
        cursor_protein.execute("SELECT * FROM {} WHERE protein_id ='{}'".format('seq_protein', seq_id))
        gene = cursor_gene.fetchall()
        protein = cursor_protein.fetchall()
        gene_seq = gene_len = pro_seq = pro_len = ''
        if gene:
            gene_seq = gene[0][1]
            gene_len = gene[0][2]
        if protein:
            pro_seq = protein[0][1]
            pro_len = protein[0][2]
        cursor_gene.close()
        cursor_protein.close()
        return gene_seq, gene_len, pro_seq, pro_len

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/prok_rna/query_seq "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_workflow",
            seq_id="YE4208",
            # seq_db_path="/mnt/ilustre/users/sanger-dev/workspace/20180831/Single_Srna_8837_fyt/Srna/output/rockhopper",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
