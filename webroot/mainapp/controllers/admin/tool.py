# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .admin import Admin
import os
import re
from biocluster.core.function import friendly_size


class ToolAction(Admin):
    def __init__(self):
        super(ToolAction, self).__init__()
        self.nav_index = 7
        self.title = "Tool工具"
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
            for d in self.list_model.get_data(self.current_page):
                if d.work_dir:
                    setattr(d, "work_dir_url", self.work_dir_url(d.work_dir))
                else:
                    setattr(d, "work_dir_url", "#")
                setattr(d, "cpu_used", d.max_cpu_use / 100)
                setattr(d, "mem_used", friendly_size(d.max_vms))
                setattr(d, "avg_cpu_used", d.average_cpu_use / 100)
                setattr(d, "avg_mem_used", friendly_size(d.average_vms))
                short_name = d.run_id.replace(d.parent_run_id, "")
                setattr(d, "short_name", short_name)
                d.run_host = re.sub(r"\.local$", "", d.run_host)
                alert_str = ""
                if d.run_spend_time > 180:
                    if abs(d.max_cpu_use / 100 - int(d.request_cpu)) > 2:
                        alert_str += "CPU申请值与实际使用值不匹配！"
                    if d.request_memory:
                        request_mem = self.memory(d.request_memory)
                        if request_mem > d.max_vms > 0:
                            if abs(d.max_vms - request_mem) / d.max_vms > 1 and abs(
                                    d.max_vms - request_mem) > 10 * 1024 * 1024 * 1024:
                                alert_str += "申请内存值与实际使用值相差过大！"
                        else:
                            if d.max_vms > 2 * 1024 * 1024 * 1024:
                                if request_mem > 0:
                                    if abs(d.max_vms - request_mem) / request_mem > 1 and abs(
                                            d.max_vms - request_mem) > 10 * 1024 * 1024 * 1024:
                                        alert_str += "申请内存值与实际使用值相差过大！"
                                else:
                                    alert_str += "请申请内存！"
                    else:
                        if d.max_vms > 2 * 1024 * 1024 * 1024:
                            alert_str += "请设置申请内存值！"
                setattr(d, "resource_alert", alert_str)
                self.list_data.append(d)
        return self.render.admin.tool(self)

    def work_dir_url(self, work_dir):
        return os.path.relpath(work_dir, self.config.WORK_DIR)

    def memory(self, value):
        """
        返回Job申请的内存大小，单位为字节数
        :return:
        """
        if value == "":
            return 0
        pattern = re.compile(r'([\d\.]+)([gmk])', re.I)
        match = re.match(pattern, value)
        if match:
            unit = match.group(2)
            if unit.upper() == "G":
                return float(match.group(1)) * 1024 * 1024 * 1024
            elif unit.upper() == "M":
                return float(match.group(1)) * 1024 * 1024
            elif unit.upper() == "K":
                return float(match.group(1)) * 1024
        else:
            pattern = re.compile(r'([\d\.]+)', re.I)
            match = re.match(pattern, value)
            if match:
                return float(match.group(1))
        return 0
