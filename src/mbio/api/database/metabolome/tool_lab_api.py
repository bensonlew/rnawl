# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from bson import SON
import datetime
import os
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from biocluster.config import Config


class ToolLabApi(Base):
    def __init__(self, bind_object):
        super(ToolLabApi, self).__init__(bind_object)
        self._project_type = 'metabolome'

    def add_circular_heatmap(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        # self.update_db_record('sg_tool_lab_circular_heatmap', main_id, status="end", main_id=main_id)
        main_collection = self.db['sg_tool_lab_circular_heatmap']
        main_collection.update({"_id": main_id}, {"$set": {'status': "end", 'main_id': main_id}})

