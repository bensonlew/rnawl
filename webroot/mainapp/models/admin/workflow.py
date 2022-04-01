# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .base import BaseModel


class WorkflowModel(BaseModel):
    def __init__(self):
        super(WorkflowModel, self).__init__()

    def get_by_workflow_id(self, workflow_id, what="*"):
        where_str = "workflow_id = '%s'" % workflow_id
        data = self.db.select(self.table, what=what, where=where_str)
        if len(data) > 0:
            return data[0]
        else:
            return None

    def get_by_run_id(self, run_id, what="*"):
        return self.get_by_workflow_id(run_id, what)

    def get_batch_list(self, batch_id, what="*"):
        where_str = "batch_id = '%s'" % batch_id
        data = self.db.select(self.table, what=what, where=where_str)
        if len(data) > 0:
            return data
        else:
            return None

    def get_report_list(self, workflow_id, end_time=None, what="*", ):
        if end_time:
            where_str = " add_time > '%s' and " % end_time
        else:
            where_str = ""
        where_str += "is_report=1 and substring(workflow_id,1,LOCATE('_',workflow_id," \
                     "POSITION('_' in workflow_id)+1)-1) = '{}' ".format(workflow_id)

        data = self.db.select(self.table, what=what, where=where_str)
        if len(data) > 0:
            return data
        else:
            return None

