# !/usr/bin/python
# -*- coding: utf-8 -*-
import types
import re
import os
import json
import math
from collections import OrderedDict
import unittest
import datetime
import glob
from bson.objectid import ObjectId
from api_base import ApiBase
import pandas as pd
import sys


class Logic(ApiBase):
    def __init__(self, bind_object):
        super(Logic, self).__init__(bind_object)


    def add_logic_base(self, logic_file):
        """
        :version:version info
        """

        df = pd.read_table(logic_file)
        df = df.fillna("")
        self.create_db_table('sgdb_logic', df.to_dict('records'))


    def get_logic_base(self, **kwargs):
        """
        :version:version info
        """

        results = self.get_tables_by_main_record("sgdb_logic", **kwargs)
        return results


    def update_logic_base(self, field_list):
        self.create_db_table('sgdb_logic', field_list)

    def merge_ensemble_uniprot(self, ensemble_file, uniprot_base_file):
        """
        合并ensemble, uniprot 表结构
        """
        pass

    def get_logic_table(self, output="logic.xls"):
        results = self.get_tables_by_main_record("sgdb_logic")
        fields =  ["field", "displayName", "ENSEMBL", "NCBI", "UNIPROT", "uniprot_use", "url", "selected", "uniprot_match"]
        with open(output, 'w') as f:
            f.write("\t".join(fields) + "\n")
            for record in results:
                values = [str(record.get(x, "")) for x in fields]
                f.write("\t".join(values) + "\n")

    def get_field2url(self):
        results = self.get_tables_by_main_record("sgdb_logic")
        field2url = {dic["field"]: dic["url"] for dic in results}
        results = self.get_tables_by_main_record("sgdb_logic")
        field2name = {dic["field"]: dic["displayName"] for dic in results}
        return field2url, field2name




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def __init__(self, method_name, params_dict=None):
        self.params_dict = params_dict
        super(TestFunction, self).__init__(methodName=method_name)
        self.toolbox = Logic(None)

    def add_logic_base(self):
        self.toolbox.add_logic_base(**self.params_dict)
    def get_logic_table(self):
        self.toolbox.get_logic_table(**self.params_dict)


if __name__ == '__main__':
    if sys.argv[1] in ["-h", "-help", "--h", "--help"]:
        print "\n".join(["add_logic_base", "get_logic_table"])
        if len(sys.argv) == 3:
            # l = Logic(None)
            help(getattr(Logic, sys.argv[2]))

    elif len(sys.argv) >= 3:
        suite = unittest.TestSuite()
        params_dict = dict()
        for par in sys.argv[2:]:
            params_dict.update({par.split("=")[0]: "=".join(par.split("=")[1:])})

        suite.addTest(TestFunction(sys.argv[1], params_dict))
        unittest.TextTestRunner(verbosity=2).run(suite)


