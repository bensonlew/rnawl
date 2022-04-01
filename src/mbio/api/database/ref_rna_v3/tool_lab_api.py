# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from bson import SON
import datetime
import os
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from bson.objectid import ObjectId
from biocluster.config import Config

class ToolLabApi(ApiBase):
    def __init__(self, bind_object):
        super(ToolLabApi, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'

    def add_diff_plot(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_diff_plot', main_id, status="end", main_id=main_id)

    def add_circularbar(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_circularbar', main_id, status="end", main_id=main_id)

    def add_vcf_distribution(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_vcf_distribution', main_id, status="end", main_id=main_id)