# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

import pandas as pd
from biocluster.core.exceptions import OptionError

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class SoftwareDatabase(ApiBase):
    """
    Api Part:
    set_db for the info of softwares and databases used by prok_rna.
    """

    def __init__(self, bind_object):
        super(SoftwareDatabase, self).__init__(bind_object)

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
            self.create_db_table("software_database", [main_info])


class TestFunction(unittest.TestCase):
    """
    This is test for the api. Just run this script to do test.
    """

    def test_mongo(test):
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        data = {
            "id": "mirna_software_database_{}".format(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "project_sn": "whole_transcriptome",
            "type": "workflow",
            "name": "whole_transcriptome.whole_transcriptome_test_api",
            "options": {},
        }
        wsheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("whole_transcriptome.software_database")
        # info_xls = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/ref_denovo/soft.xls"
        info_xls = "/mnt/lustre/users/sanger/sg-users/qinjincheng/whole_transcriptome/soft.xls"
        wf.test_api.add_swdb_info(info_xls)


if __name__ == "__main__":
    unittest.main()
