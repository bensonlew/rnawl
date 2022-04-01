# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .base import BaseModel


class ToolModel(BaseModel):

    def __init__(self):
        super(ToolModel, self).__init__()

    def get_workflow_tools(self, workflow_id, success=True):
        where_str = "run_id like '%s.%%'" % workflow_id
        if success:
            where_str += " and success=1"
        data = self.db.select(self.table, where=where_str)
        if len(data) > 0:
            return data
        else:
            return None

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
