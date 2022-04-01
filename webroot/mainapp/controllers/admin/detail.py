# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .admin import Admin
from biocluster.core.function import get_clsname_form_path, friendly_size
import os
import re


class DetailAction(Admin):
    def __init__(self):
        super(DetailAction, self).__init__()
        self.check_required(["run_id"])
        self.nav_index = 0
        self.title = "运行详情"
        self.type = "workflow"
        self.data_index = 0
        self.data_list = []
        self.show_type = "sub"
        self.is_batch = False
        self.start = 0
        self.end = 0

    def GET(self):
        if self.params.show == "report":
            self.show_type = "report"
            data = self.get_all_childs(get_sub=False)
            if data:
                setattr(data, "data_type", "workflow")
                if self.type == "workflow":
                    if data.is_end == 1:
                        batch_list = self.get_report(self.params.run_id, data.end_time)
                    else:
                        batch_list = None
                    if batch_list:
                        setattr(data, "reports", batch_list)
                self.get_data([data], 0)
        elif self.params.show == "batch":
            self.show_type = "batch"
            data = self.get_all_childs(get_sub=False)
            if data:
                setattr(data, "data_type", "workflow")
                if self.type == "workflow" and self.is_batch:
                    batch_list = self.get_batch(self.params.run_id)
                    if batch_list:
                        setattr(data, "batchs", batch_list)
                self.get_data([data], 0)

        else:
            self.show_type = "sub"
            data = self.get_all_childs()
            if data:
                self.get_data([data], 0)
        return self.render.admin.detail(self)

    @property
    def workflow_id(self):
        if self.params.run_id:
            return self.params.run_id.split(".")[0]
        else:
            return None

    def work_dir(self, data):
        if data.work_dir:
            return data.work_dir
        else:
            if data.type and data.path and data.run_time:
                timestr = data.run_time.strftime('%Y%m%d')
                work_dir = self.config.WORK_DIR + "/" + timestr + "/" + self.__get_min_name(data.type, data.path) + "_" + data.workflow_id
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
            if re.search(b + "$", class_name):
                return re.sub((b + "$"), "", class_name)
        return class_name

    def work_dir_url(self, work_dir):
        return os.path.relpath(work_dir, self.config.WORK_DIR)

    def get_batch(self, run_id):
        model = self.get_model("workflow")
        data = model.get_batch_list(run_id, what="id,client,workflow_id,add_time,is_end,end_time,"
                                                 "server,is_error,error,pid,paused,instant,has_run,"
                                                 "run_time,path,type,batch_id,cpu_used,memory_used,"
                                                 "batch,work_dir,rerun,rerun_time,cluster")
        data_list = []
        if data:
            for d in data:
                setattr(d, "data_type", "workflow")
                work_dir = self.work_dir(d)
                if work_dir != "#":
                    setattr(d, "work_dir", work_dir)
                else:
                    setattr(d, "work_dir", "#")
                    setattr(d, "work_dir_url", "#")
                data_list.append(d)
        return data_list

    def get_report(self, run_id, end_time):
        model = self.get_model("workflow")
        data = model.get_report_list(run_id, end_time, what="id,client,workflow_id,add_time,is_end,end_time,"
                                                            "server,is_error,error,pid,paused,instant,has_run,"
                                                            "run_time,path,type,batch_id,cpu_used,memory_used,"
                                                            "batch,work_dir,rerun,rerun_time,cluster")
        data_list = []
        if data:
            for d in data:
                setattr(d, "data_type", "workflow")
                work_dir = self.work_dir(d)
                if work_dir != "#":
                    setattr(d, "work_dir", work_dir)
                else:
                    setattr(d, "work_dir", "#")
                    setattr(d, "work_dir_url", "#")
                data_list.append(d)
        return data_list

    def get_data(self, data=None, index=0):
        for d in data:
            self.data_index += 1
            setattr(d, "i_index", self.data_index)
            if index > 0:
                setattr(d, "parent_index", index)
            m = None
            t = None
            c = None
            b = None
            r = None
            if hasattr(d, "modules"):
                m = d.modules
                delattr(d, "modules")
            if hasattr(d, "tools"):
                t = d.tools
                delattr(d, "tools")
            if hasattr(d, "commands"):
                c = d.commands
                delattr(d, "commands")
            if hasattr(d, "batchs"):
                b = d.batchs
                delattr(d, "batchs")
            if hasattr(d, "reports"):
                r = d.reports
                delattr(d, "reports")
            if hasattr(d, "work_dir"):
                setattr(d, "work_dir_url", self.work_dir_url(d.work_dir))
            if d.data_type == "workflow":
                if d.run_time:
                    setattr(d, "queue_spend_time", (d.run_time - d.add_time).seconds)
                else:
                    setattr(d, "queue_spend_time", 0)
                if d.end_time:
                    setattr(d, "run_spend_time", (d.end_time - d.run_time).seconds)
                else:
                    setattr(d, "run_spend_time", 0)
            if d.data_type == "module":
                setattr(d, "queue_spend_time", "-")
                if d.end_time:
                    setattr(d, "run_spend_time", (d.end_time - d.start_time).seconds)
                else:
                    setattr(d, "run_spend_time", 0)
            if d.data_type == "command":
                setattr(d, "queue_spend_time", "-")
                if d.end_time:
                    setattr(d, "run_spend_time", (d.end_time - d.start_time).seconds)
                else:
                    setattr(d, "run_spend_time", 0)
                setattr(d, "cpu_used", d.max_cpu_use/100)
                setattr(d, "mem_used", friendly_size(d.max_vms))
                setattr(d, "avg_cpu_used", d.average_cpu_use / 100)
                setattr(d, "avg_mem_used", friendly_size(d.average_vms))
            if d.data_type == "tool":
                setattr(d, "cpu_used", d.max_cpu_use/100)
                setattr(d, "mem_used", friendly_size(d.max_vms))
                setattr(d, "avg_cpu_used", d.average_cpu_use / 100)
                setattr(d, "avg_mem_used", friendly_size(d.average_vms))
                d.run_host = re.sub(r"\.local$", "", d.run_host)
                alert_str = ""
                if d.run_spend_time > 180:
                    if abs(d.max_cpu_use/100 - int(d.request_cpu)) > 2:
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
                                    if abs(d.max_vms - request_mem) / request_mem > 1 and abs(d.max_vms - request_mem) > 10 * 1024 * 1024 * 1024:
                                        alert_str += "申请内存值与实际使用值相差过大！"
                                else:
                                    alert_str += "请申请内存！"
                    else:
                        if d.max_vms > 2 * 1024 * 1024 * 1024:
                            alert_str += "请设置申请内存值！"
                setattr(d, "resource_alert", alert_str)

            self.data_list.append(d)

            if m:
                self.get_data(m, d.i_index)
            if t:
                self.get_data(t, d.i_index)
            if c:
                self.get_data(c, d.i_index)
            if b:
                self.get_data(b, d.i_index)
            if r:
                self.get_data(r, d.i_index)

    def get_all_childs(self, get_sub=True):
        model = self.get_model("workflow")
        data = model.get_by_run_id(self.params.run_id, what="id,client,workflow_id,add_time,is_end,end_time,"
                                                            "server,is_error,error,pid,paused,instant,has_run,"
                                                            "run_time,path,type,batch_id,cpu_used,memory_used,"
                                                            "batch,work_dir,rerun,rerun_time,cluster")
        if data:
            setattr(data, "data_type", "workflow")
            work_dir = self.work_dir(data)
            if work_dir != "#":
                setattr(data, "work_dir", work_dir)
            else:
                setattr(data, "work_dir", "#")
                setattr(data, "work_dir_url", "#")
            self.type = "workflow"
            if data.batch == 1:
                self.is_batch = True
            if get_sub:
                data1 = self.get_child_modules(data.workflow_id, data.cluster)

                if data1:
                    setattr(data, "modules", data1)
                data2 = self.get_child_tools(data.workflow_id, data.cluster)
                if data2:
                    setattr(data, "tools", data2)
        else:
            d = model.get_by_run_id(self.workflow_id, what="cluster")
            cluster = d.cluster
            model = self.get_model("module")
            data = model.get_by_run_id(self.params.run_id)
            if data:
                setattr(data, "data_type", "module")
                self.type = "module"
                setattr(data, "short_name", data.run_id)
                setattr(data, "cluster", cluster)
                if get_sub:
                    data1 = self.get_child_modules(data.run_id, cluster)
                    if data1:
                        setattr(data, "modules", data1)
                    data2 = self.get_child_tools(data.run_id, cluster)
                    if data2:
                        setattr(data, "tools", data2)
            else:
                model = self.get_model("tool")
                data = model.get_by_run_id(self.params.run_id)
                if data:
                    setattr(data, "data_type", "tool")
                    self.type = "tool"
                    setattr(data, "short_name", data.run_id)
                    setattr(data, "cluster", cluster)
                    if get_sub:
                        data1 = self.get_child_commands(data.run_id, cluster, data.work_dir)
                        if data1:
                            setattr(data, "tools", data1)

        return data

    def get_child_modules(self, run_id, cluster):
        model = self.get_model("module")
        data = model.get_by_parent_id(run_id)
        data_list = []
        if data:
            for d in data:
                setattr(d, "data_type", "module")
                short_name = d.run_id.replace(run_id, "")
                setattr(d, "short_name", short_name)
                setattr(d, "cluster", cluster)
                # data1 = self.get_child_modules(d.run_id)
                # if data1:
                #     setattr(d, "modules", data1)
                # data2 = self.get_child_tools(d.run_id)
                # if data2:
                #     setattr(d, "tools", data2)
                data_list.append(d)
        return data_list

    def get_child_tools(self, run_id, cluster):
        model = self.get_model("tool")
        data = model.get_by_parent_id(run_id)
        data_list = []
        if data:
            for d in data:
                setattr(d, "data_type", "tool")
                short_name = d.run_id.replace(run_id, "")
                setattr(d, "short_name", short_name)
                setattr(d, "cluster", cluster)
                data1 = self.get_child_commands(d.run_id, cluster, d.work_dir)
                if data1:
                    setattr(d, "commands", data1)
                data_list.append(d)
        return data_list

    def get_child_commands(self, run_id, cluster, parent_dir=None):
        model = self.get_model("command")
        data = model.get_by_parent_id(run_id)
        data_list = []
        if data:
            for d in data:
                setattr(d, "data_type", "command")
                setattr(d, "cluster", cluster)
                if parent_dir:
                    setattr(d, "work_dir", parent_dir)
                data_list.append(d)
        return data_list

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
            unit = str(match.group(2))
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
