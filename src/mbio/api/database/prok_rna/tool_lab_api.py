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
        self._project_type = 'prok_rna'

    def add_upset_venn(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_upset_venn', main_id, status="end", main_id=main_id)

    def add_primer(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_seqdown', main_id, status="end", main_id=main_id)

    def add_network(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_network', main_id, status="end", main_id=main_id)

    def add_geneset_ppi(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_geneset_ppi', main_id, status="end", main_id=main_id)

    def add_go_circ(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_go_circ', main_id, status="end", main_id=main_id)

    def add_wgcna(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_wgcna', main_id, status="end", main_id=main_id)