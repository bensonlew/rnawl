# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .base import BaseModel


class ModuleModel(BaseModel):
    def __init__(self):
        super(ModuleModel, self).__init__()

    def get_by_run_id(self, run_id, what="*"):
        where_str = "run_id = '%s'" % run_id
        data = self.db.select(self.table, what=what, where=where_str)
        if len(data) > 0:
            return data[0]
        else:
            return None

    def get_by_parent_id(self, parent_run_id, what="*"):
        where_str = "parent_run_id = '%s'" % parent_run_id
        data = self.db.select(self.table, what=what, where=where_str)
        if len(data) > 0:
            return data
        else:
            return None





