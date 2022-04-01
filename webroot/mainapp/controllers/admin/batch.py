# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .workflow import WorkflowAction


class BatchAction(WorkflowAction):
    def __init__(self):
        super(BatchAction, self).__init__()
        self.nav_index = 5
        self.title = "批处理流程"
        self.list_data = []
        self._table = "workflow"

    def GET(self):
        self.list_model.pagesize = 50
        self.list_model.set_conditions(where="batch = 1", what="id,client,workflow_id,add_time,is_end,end_time,server,"
                                                               "is_error,error,pid,paused,instant,has_run,run_time,"
                                                               "path,type,batch_id,cpu_used,memory_used,batch,work_dir,"
                                                               "cluster")

        if self.list_model.total > 0:
            for data in self.list_model.get_data(self.current_page):
                work_dir = self.work_dir(data)
                if work_dir != "#":
                    setattr(data, "work_dir", work_dir)
                    setattr(data, "work_dir_url", self.work_dir_url(work_dir))
                else:
                    setattr(data, "work_dir", "未知")
                    setattr(data, "work_dir_url", "#")
                self.list_data.append(data)
        return self.render.admin.batch(self)


