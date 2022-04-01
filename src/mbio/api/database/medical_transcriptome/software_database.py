#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.medical_transcriptome.api_base import ApiBase
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
        self._project_type = "medical_transcriptome"

    def add_swdb_info(self, info_xls):
        if not os.path.isfile(info_xls):
            self.bind_object.set_error("Error: the value passed of info_xls is not file!", code="53703502")
        df = pd.read_table(info_xls)
        for i in df.index:
            line_dict = df.loc[i]
            main_info = {
                "software_database": line_dict["software_database"],
                "source": line_dict["source"],
                "usage": line_dict["usage"],
                "version": line_dict["version"],
                "task_version" : "v1",
                # "task_version" : line_dict["task_version"],
            }
            self.create_db_table("sg_software_database", [main_info])

    def update_record_by_dict(self, collection_name, query_dict, insert_dict):
        conn = self.db[collection_name]
        conn.update(query_dict, {"$set": insert_dict}, multi=True, upsert=True)

    def download_db(self, collection_name):
        conn = self.db[collection_name]
        records = conn.find({})

        print "["
        print ",\n".join([str(record) for record in records])
        print "]"

    def update_swdb_info(self, task_version=None):
        query_dict = {}
        insert_dict = {"task_version": task_version}
        self.update_record_by_dict("sg_software_database", query_dict, insert_dict)

class TestFunction(unittest.TestCase):
    """
    This is test for the api. Just run this script to do test.
    """
    def test_mongo(test):
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        data = {
            "id": "mirna_software_database_{}".format(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "project_sn": "medical_transcriptome",
            "type": "workflow",
            "name": "medical_transcriptome.medical_transcriptome_test_api",
            "options": {},
        }
        wsheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("medical_transcriptome.software_database")
        info_xls = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/ref_rna_v3/medical_transcriptome_new.xls"
        wf.test_api.add_swdb_info(info_xls)
        # wf.test_api.download_db("sg_software_database")

if __name__ == "__main__":
    unittest.main()
