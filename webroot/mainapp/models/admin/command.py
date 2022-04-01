# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .base import BaseModel


class CommandModel(BaseModel):
    def __init__(self):
        super(CommandModel, self).__init__()

    def get_by_parent_id(self, parent_run_id, what="*"):
        where_str = "parent_run_id = '%s'" % parent_run_id
        data = self.db.select(self.table, what=what, where=where_str)
        if len(data) > 0:
            return data
        else:
            return None
