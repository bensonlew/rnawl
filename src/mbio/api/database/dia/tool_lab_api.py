# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from bson import SON
import datetime
import os
from mbio.api.database.dia.api_base import ApiBase
from bson.objectid import ObjectId
from biocluster.config import Config

class ToolLabApi(ApiBase):
    def __init__(self, bind_object):
        super(ToolLabApi, self).__init__(bind_object)
        self._project_type = 'dia'

    def add_upset_venn(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_upset_venn', main_id, status="end", main_id=main_id)

    def add_diff_plot(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_diff_plot', main_id, status="end", main_id=main_id)

    def add_diff_ma(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_diff_ma', main_id, status="end", main_id=main_id)

    def add_diff_radar(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_diff_radar', main_id, status="end", main_id=main_id)

    def add_pathway_network(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_pathway_network', main_id, status="end", main_id=main_id)

    def add_circularbar(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_circularbar', main_id, status="end", main_id=main_id)

    def add_circular_heatmap(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        self.update_db_record('sg_tool_lab_circular_heatmap', main_id, status="end", main_id=main_id)

