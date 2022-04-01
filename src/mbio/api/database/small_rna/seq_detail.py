# !/usr/bin/python
# -*- coding: utf-8 -*-
import re
import os
import pandas as pd
import unittest
from mbio.api.database.small_rna.api_base import ApiBase
import sqlite3


class SeqDetail(ApiBase):
    def __init__(self, bind_object):
        super(SeqDetail, self).__init__(bind_object)

    @staticmethod
    def parse_seq_file(seq_file):
        """
        generator for parsing sequence file
        :param seq_file: fasta sequence file
        :param match: regexp pattern for parsing line startswith '>'
        :return: seq_tag, sequence
        """
        with open(seq_file, 'r') as f:
            j = 0
            seq_tag = tuple()
            sequence = ""
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        yield seq_tag, sequence
                    seq_tag = line.strip()
                    sequence = ''
                else:
                    sequence += line.strip()
            else:
                yield seq_tag, sequence

    def build_seq_database(self, db_path, cds_file, pep_file, transcript, trans2unigene, task_id='denovo_rna_v2'):
        """
        build sequence db
        :param db_path: abs path of seq db to be build.
        :param cds_file:
        :param pep_file:
        :param transcript: transcript fasta from trinity
        :param trans2unigene: file with two columns [seq_id, 'yes/No']
        :param task_id: task_id value, db_path will be saved in sg_task
        :return:
        """
        df = pd.read_table(trans2unigene, sep='\t', header=None)
        t2g_dict = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
        is_unigene_trans = list(df[df.iloc[:, 2] == 'yes'].iloc[:, 0])
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        # Create cds info table
        table_name = "seq_cds"
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} (t_id text, g_id text, cds_pos text, cds text, pep text)'.format(table_name))
        # add cds and pep detail
        # save pep info into dictionary first
        pep_parser = self.parse_seq_file(pep_file)
        pep_dict = dict()
        for pep_desc, pep_seq in pep_parser:
            pep_dict[pep_desc] = pep_seq
        match = re.compile(r'>(.*?)::(.*?)::.*type:(.*?)\s+len:\d+.*:(\d+-\d+\(.\)).*').match
        cds_parser = self.parse_seq_file(cds_file)
        for cds_desc, cds_seq in cds_parser:
            # ('TRINITY_DN1001_c0_g1', 'TRINITY_DN1001_c0_g1_i1', 'complete', '92-931(+)')
            g_id, t_id, _type, cds_pos = match(cds_desc).groups()
            pep_seq = pep_dict[cds_desc]
            if t_id not in is_unigene_trans:
                g_id += ": not_unigene"
            cursor.execute("INSERT INTO {} VALUES ('{}','{}', '{}', '{}', '{}')".format(
                "seq_cds", t_id, g_id, cds_pos, cds_seq, pep_seq))
        # add sequence detail
        # Create sequence table for gene or transcript
        table_name = "transcript_gene"
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} (t_id text, g_id text, sequence text)'.format(table_name))
        parser = self.parse_seq_file(transcript)
        match_id = re.compile(r'>([^\s]+)').match
        for seq_tag, seq in parser:
            # ('TRINITY_DN1001_c0_g1', 'TRINITY_DN1001_c0_g1_i1', 'complete', '92-931(+)')
            t_id = match_id(seq_tag).groups()[0]
            g_id = t2g_dict[t_id]
            if t_id not in is_unigene_trans:
                g_id += ": not_unigene"
            cursor.execute("INSERT INTO {} VALUES ('{}','{}', '{}')".format(
                "transcript_gene", t_id, g_id, seq))
        #
        conn.commit()
        conn.close()
        # finish local database building for sequence
        #self.update_record_by_task_id('sg_task', task_id, seq_db=db_path)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/api/database/' \
               'denovo_rna_v2/test_files'
    toolbox = SeqDetail(None)

    def test_build_seq_databse(self):
        cds = os.path.join(self.test_dir, 'Trinity.fasta.transdecoder.cds')
        pep = os.path.join(self.test_dir, 'Trinity.fasta.transdecoder.pep')
        fasta = os.path.join(self.test_dir, 'Trinity.fasta')
        trans2unigene = os.path.join(self.test_dir, 'Trinity.fasta_t2g2u')
        seq_db = '/mnt/ilustre/users/sanger-dev/sg-users/deqing/seqs.db'
        self.toolbox.build_seq_database(seq_db, cds, pep, fasta, trans2unigene)
        conn = sqlite3.connect(seq_db)
        cursor = conn.cursor()
        table_name = "seq_cds"
        query_id = "TRINITY_DN1001_c0_g1"
        query_id = "TRINITY_DN1004_c0_g1"
        cursor.execute("SELECT * FROM {} WHERE g_id='{}'".format(table_name, query_id))
        results = cursor.fetchall()
        # print(results)
        if results:
            data = list()
            for each in results:
                print(each)
        else:
            print("{} not found".format(query_id))
        #
        table_name = 'transcript_gene'
        cursor.execute("SELECT sequence FROM {} WHERE g_id='{}'".format(table_name, query_id))
        print(cursor.fetchall()[0][0])
        cursor.close()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_build_seq_databse'))
    unittest.TextTestRunner(verbosity=2).run(suite)

