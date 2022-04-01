#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.prok_rna.api_base import ApiBase
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest

class SoftwareDatabase(ApiBase):
    """
    Api Part:
    set_db for the info of softwares and databases used by prok_rna.
    """
    def __init__(self, bind_object):
        super(SoftwareDatabase, self).__init__(bind_object)
        self._project_type = "small_rna"

    def add_swdb_info(self, info_xls):
        if not os.path.isfile(info_xls):
            raise OptionError("Error: the value passed of info_xls is not file!")
        df = pd.read_table(info_xls)
        for i in df.index:
            line_dict = df.loc[i]
            main_info = {
                "software_database": line_dict["software_database"],
                "source": line_dict["source"],
                "usage": line_dict["usage"],
                "version": line_dict["version"],
            }
            self.create_db_table("sg_software_database", [main_info])

class TestFunction(unittest.TestCase):
    """
    This is test for the api. Just run this script to do test.
    """
    def test_mongo(test):
        from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        data = {
            "id": "mirna_software_database_{}".format(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "project_sn": "small_rna",
            "type": "workflow",
            "name": "small_rna.smallrna_test_api",
            "options": {},
        }
        wsheet = Sheet(data=data)
        wf = SmallRnaTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("small_rna.software_database")
        info_xls = "/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/software_database.txt"
        wf.test_api.add_swdb_info(info_xls)

if __name__ == "__main__":
    unittest.main()
