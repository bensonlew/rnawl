# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .admin import Admin
from biocluster.core.function import get_clsname_form_path, pretty_date, friendly_time, get_error_str
import re
import os
import datetime
import time

class WorkflowAction(Admin):
    
    def __init__(self):
        super(WorkflowAction, self).__init__()
        self.nav_index = 4
        self.title = "Workflow流程"
        self.list_data = []

    def GET(self):
        self.list_model.pagesize = 50
        where_str = None

        if self.params.data_type:
            if self.params.table_search:
                if self.params.data_type == 'run_id':
                    where_str = "workflow_id='%s'" % self.params.table_search
                if self.params.data_type == 'path':
                    where_str = "path='%s'" % self.params.table_search
            elif self.params.time_from or self.params.time_to:
                if self.params.data_type in ["add_time", "end_time", "run_time"]:
                    key = self.params.data_type
                else:
                    key = "add_time"
                if self.params.time_from:
                    where_str = "%s>='%s' " % (key, self.params.time_from)
                if self.params.time_to:
                    if where_str:
                        where_str += "and %s<='%s' " % (key, self.params.time_to)
                    else:
                        where_str = "%s<='%s' " % (key, self.params.time_to)
        if self.params.status:
            if self.params.status == 'running':
                where_str = "has_run = 1 and is_end = 0 and paused = 0"
            elif self.params.status == 'queue':
                where_str = "has_run = 0"
            elif self.params.status == 'error':
                where_str = "has_run = 1 and is_error = 1 and is_end = 1"
            elif self.params.status == 'success':
                where_str = "has_run = 1 and is_error = 0 and is_end = 1"
            elif self.params.status == 'paused':
                where_str = "has_run = 1 and paused = 1 and is_end = 0"
        self.list_model.set_conditions(where=where_str, what="id,client,workflow_id,add_time,is_end,end_time,server,"
                                                             "is_error,error,pid,paused,instant,has_run,run_time,"
                                                             "path,type,batch_id,cpu_used,memory_used,batch,work_dir,"
                                                             "rerun, rerun_time, cluster")
        if self.list_model.total > 0:
            for data in self.list_model.get_data(self.current_page):
                work_dir = self.work_dir(data)
                if work_dir != "#":
                    setattr(data, "work_dir", work_dir)
                    setattr(data, "work_dir_url", self.work_dir_url(work_dir))
                else:
                    setattr(data, "work_dir", "未知")
                    setattr(data, "work_dir_url", "#")
                if data.cpu_used is None:
                    data.cpu_used = 0
                if data.memory_used is None:
                    data.memory_used = 0
                now = datetime.datetime.now()
                if data.rerun == 1 and data.rerun_time:
                    start = data.rerun_time
                else:
                    start = data.add_time
                if data.is_error == 1:
                    data.error = get_error_str(data.error)
                if data.is_end == 1 or data.is_error == 1:
                    time_spend = int(time.mktime(data.end_time.timetuple())) - \
                                 int(time.mktime(data.add_time.timetuple()))
                    data.end_time_str = pretty_date(now - data.end_time)
                else:
                    if now > start:
                        time_spend = int(time.mktime(now.timetuple())) - \
                                 int(time.mktime(start.timetuple()))
                    else:
                        time_spend = 0
                    data.end_time_str = " "
                data.spend_time = friendly_time(time_spend)
                data.add_time = pretty_date(now-start)
                self.list_data.append(data)
        return self.render.admin.workflow(self)

    def work_dir(self, data):
        if data.work_dir:
            return data.work_dir
        else:
            if data.type and data.path and data.run_time:
                timestr = data.run_time.strftime('%Y%m%d')
                work_dir = self.config.WORK_DIR + "/" + timestr + "/" + self.__get_min_name(data.type, data.path) \
                           + "_" + data.workflow_id
                return work_dir
            else:
                return "#"

    @staticmethod
    def __get_min_name(typename, name):
        if typename != "workflow":
            return "Singal"
        class_name = get_clsname_form_path(name)
        base = ["Basic", "Module", "Tool", "Agent", "Workflow"]
        for b in base:
            if re.search(b+"$", class_name):
                return re.sub((b + "$"), "", class_name)
        return class_name

    def work_dir_url(self, work_dir):
        return os.path.relpath(work_dir, self.config.WORK_DIR)
