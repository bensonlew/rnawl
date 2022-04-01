# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .admin import Admin
import os


class ModuleAction(Admin):
    def __init__(self):
        super(ModuleAction, self).__init__()
        self.nav_index = 6
        self.title = "Module模块"
        self.list_data = []

    def GET(self):
        self.list_model.pagesize = 50
        where_str = None

        if self.params.data_type:
            if self.params.table_search:
                if self.params.data_type == 'run_id':
                    where_str = "run_id='%s'" % self.params.table_search
                if self.params.data_type == 'path':
                    where_str = "path='%s'" % self.params.table_search
            elif self.params.time_from or self.params.time_to:
                if self.params.data_type == "end_time":
                    key = self.params.data_type
                else:
                    key = "start_time"
                if self.params.time_from:
                    where_str = "%s>='%s' " % (key, self.params.time_from)
                if self.params.time_to:
                    if where_str:
                        where_str += "and %s<='%s' " % (key, self.params.time_to)
                    else:
                        where_str = "%s<='%s' " % (key, self.params.time_to)
        self.list_model.set_conditions(where=where_str)
        if self.list_model.total > 0:
            for data in self.list_model.get_data(self.current_page):
                if data.work_dir:
                    setattr(data, "work_dir_url", self.work_dir_url(data.work_dir))
                else:
                    setattr(data, "work_dir_url", "#")
                setattr(data, "workflow_id", self.workflow_id(data.run_id))
                self.list_data.append(data)
        return self.render.admin.module(self)

    def work_dir_url(self, work_dir):
        return os.path.relpath(work_dir, self.config.WORK_DIR)

    @staticmethod
    def workflow_id(run_id):
        if run_id:
            return run_id.split(".")[0]
        else:
            return None

