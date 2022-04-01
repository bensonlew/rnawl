# !/usr/bin/python
# -*- coding: utf-8 -*-
from bson.objectid import ObjectId
import types
import re
import os
import json
import pandas as pd
import numpy as np
from scipy import stats
import math
from sklearn import decomposition
import fastcluster as hclust
from collections import OrderedDict
import unittest
import datetime
from bson.son import SON
import gevent
import glob
from biocluster.api.database.base import Base, report_check
from api_base import ApiBase
import sqlite3


class ProteinSeq(ApiBase):
    def __init__(self, bind_object):
        super(ProteinSeq, self).__init__(bind_object)

    @staticmethod
    def parse_seq_file(seq_file, match=re.compile(r'>([^\s]+)').match):
        """
        generator for parsing sequence file
        :param seq_file: fasta sequence file
        :param match: regexp pattern for parsing line startswith '>'
        :return: seq_tag, sequence
        """
        with open(seq_file, 'r') as f:
            j = 0
            seq_tag = ""
            sequence = ""
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        yield seq_tag, sequence
                        #print seq_tag, sequence
                    seq_tag = match(line).groups()[0]
                    sequence = ''
                else:
                    sequence += line.strip()
            else:
                yield seq_tag, sequence

    def build_seq_database(self, db_path, pep_file, task_id='denovo_rna_v2'):
        """
        build sequence db
        """
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        # add pep detail
        table_name = "seq_protein"
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} (pep_id text, pep_seq text)'.format(table_name))
        match = re.compile(r'>(.*?)\s+.*').match
        pep_parser = self.parse_seq_file(pep_file, match)
        for pep in pep_parser:
            #print pep[0]
            pep_id = pep[0]
            pep_seq = pep[1]
            cursor.execute("INSERT INTO {} VALUES ('{}','{}')".format(
                "seq_protein",pep_id, pep_seq))
        conn.commit()
        conn.close()
        # finish local database building for sequence
        # self.update_record_by_task_id('sg_task', task_id, seq_db=db_path)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/api/database/' \
               'denovo_rna_v2/test_files'
    toolbox = ProteinSeq(None)

    def test_build_seq_databse(self):
        pep = "/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/itraq_test/DB1.fasta"
        seq_db = "/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/itraq_test/seq_db.sqlite3"
        conn = sqlite3.connect(seq_db)
        cursor = conn.cursor()
        table_name = "seq_protein"
        query_id = "P93206"
        cursor.execute("SELECT * FROM {} WHERE pep_id='{}'".format(table_name, query_id))
        results = cursor.fetchall()
        # print(results)
        if results:
            data = list()
            for each in results:
                print(each)
        else:
            print("{} not found".format(query_id))
        cursor.close()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_build_seq_databse'))
    unittest.TextTestRunner(verbosity=2).run(suite)

