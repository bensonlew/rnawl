# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .admin import Admin
from biocluster.core.function import get_clsname_form_path
import os
import re
import sys
import time
import datetime
reload(sys)
sys.setdefaultencoding('utf8')


class LogModel(object):
    def __init__(self):
        self.current_page = 1
        self.total = 1
        self.pagesize = 500

    @property
    def total_page(self):
        if self.total > 0:
            return int(self.total/self.pagesize) + 1


class LogAction(Admin):
    def __init__(self):
        super(LogAction, self).__init__()
        self.check_required(["id", "workflow_id", "run_id"], relation="or")
        self.nav_index = 4
        self.title = "日志查看"
        self.log_path = []
        self.workflow_id = None
        self.print_lines = 0
        self.id = None
        self.filer_key = None
        self.start_date = None
        self.end_date = None
        self.time = None
        self.area = None
        self.type = "workflow"
        self._is_post = False
        self._list_model = LogModel()
        self.cluster = self.params.cluster if self.params.cluster else "sanger"

    @property
    def start(self):
        return (self.current_page - 1) * self.list_model.pagesize

    @property
    def end(self):
        return self.current_page * self.list_model.pagesize

    def GET(self):
        self.list_model.current_page = self.current_page
        return self.get_file_path()

    @staticmethod
    def __get_min_name(name):
        class_name = get_clsname_form_path(name)
        base = ["Basic", "Module", "Tool", "Agent", "Workflow"]
        for b in base:
            if re.search(b + "$", class_name):
                return re.sub((b + "$"), "", class_name)
        return class_name

    def get_file_path(self):
        file_path = []
        if self.params.type == "tool":
            model = self.get_model("tool")
            self.type = "Tool"
        elif self.params.type == "command":
            model = self.get_model("command")
            self.type = "Command"
        else:
            model = self.get_model("workflow")
        if self.params.id:
            data = model.get_by_id(self.params.id)
        elif self.params.workflow_id:
            data = model.get_by_workflow_id(self.params.workflow_id)
        else:
            data = model.get_by_run_id(self.params.run_id)
        if data:
            if hasattr(data, "workflow_id"):
                self.workflow_id = data.workflow_id
            elif hasattr(data, "run_id"):
                self.workflow_id = data.run_id
            else:
                self.workflow_id = data.name
            if hasattr(data, "start_time"):
                self.start_date = data.start_time.strftime("%Y-%m-%d %H:%M:%S")
            else:
                self.start_date = data.add_time.strftime("%Y-%m-%d %H:%M:%S")
            if hasattr(data, "end_date"):
                self.end_date = data.end_time.strftime("%Y-%m-%d %H:%M:%S")
            else:
                self.end_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            self.id = data.id
            if self.params.type == "tool":
                name = self.__get_min_name(data.path)
                file_path.append(os.path.join(data.work_dir, "%s_%s.err" % (name, data.job_id)))
                file_path.append(os.path.join(data.work_dir, "%s_%s.out" % (name, data.job_id)))
                print file_path
            elif self.params.type == "command":
                name = data.name
                model1 = self.get_model("tool")
                data1 = model1.get_by_run_id(data.parent_run_id)
                file_path.append(os.path.join(data1.work_dir, "%s.o" % name))
            else:
                if data.has_run:
                    if data.rerun:
                        timestr = data.rerun_time.strftime('%Y%m%d')
                    else:
                        timestr = data.run_time.strftime('%Y%m%d')
                    log_dir = os.path.join(self.config.wpm_log_file, timestr)
                    file_path.append(os.path.join(log_dir, "%s_%s.log" %
                                                  (data.workflow_id, re.sub(r"\.local$", "", data.server))))
                else:
                    return self.failed("尚未开始运行，无法查看日志！")
            passed = False
            for f in file_path:
                if os.path.exists(f):
                    self.log_path.append(f)
                    self.list_model.total += int(os.popen('wc -l %s' % f).read().split()[0])
                    passed = True

            if not passed:
                return self.failed("日志文件不存在，无法查看！", "workflow")
            else:
                if self._is_post:
                    return self.code
                else:
                    return self.render.admin.log(self)
        else:
            return self.failed("没有找到相关数据！")

    def POST(self):
        self.list_model.current_page = self.current_page
        self._is_post = True
        if self.params.key:
            self.filer_key = self.params.key
        else:
            self.check_required(["time", "area"])
            self.time = time.mktime(time.strptime(self.params.time, '%Y-%m-%d %H:%M:%S'))
            self.area = int(self.params.area)
        return self.get_file_path()

    @property
    def code(self):
        data = ""
        in_area = False
        count = 0
        for log in self.log_path:
            with open(log, "r") as f:
                work_dir = os.environ['HOME']
                line = f.readline()
                trace_start = False
                while line:
                    if count < self.start:
                        count += 1
                        line = f.readline()
                        continue
                    elif count > self.end:
                        break
                    if self.filer_key:
                        if line.find(self.filer_key) < 0:
                            line = f.readline()
                            continue
                    m = re.match(r'^(\d{4}\-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', line)
                    if m:
                        if self.time:
                            if abs(self.time - time.mktime(time.strptime(m.group(1), '%Y-%m-%d %H:%M:%S'))) > self.area:
                                line = f.readline()
                                in_area = False
                                continue
                            else:
                                in_area = True

                        # if not self.start_date:
                        #     self.start_date = m.group(1)
                        # self.end_date = m.group(1)
                    else:
                        if self.time:
                            if not in_area:
                                line = f.readline()
                                continue
                        self.print_lines += 1
                    if trace_start and re.match(r"^\w+:", line):
                        trace_start = False
                    if re.match(r"^Traceback", line):
                        trace_start = True
                    if trace_start:
                        m = re.match(r"^\s+File \"(.+)\", line (\d+), (.*)$", line)
                        if m:
                            link_str = "<a target=\"_blank\" href=\"pymodule?cluster=%s&path=%s&line=%s\">%s</a>" % \
                                       (self.cluster, os.path.relpath(m.group(1), work_dir), m.group(2), m.group(1))
                            data += "  File \"%s\", line %s, %s\n" % (link_str, m.group(2), m.group(3))
                        else:
                            data += line
                    else:
                        m = re.match(r"^\d{4}\-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\s+\S+\((\S+)\)\s+(.*)", line)
                        if m:
                            link_str = "<a target=\"_blank\" href=\"detail?run_id=%s\">%s</a>" % (m.group(1), m.group(1))
                            data += re.sub(m.group(1), link_str, line, 1)
                            if m.group(2).find("get data to url") > -1 or m.group(2).find("post data to url"):
                                data += self.read_send(f)
                        else:
                            data += line
                    line = f.readline()
                    count += 1
        return data

    @staticmethod
    def read_send(f):
        data = ""
        while True:
            last_tell = f.tell()
            line = f.readline()
            mm = re.match(r"^\s*\n$", line)
            if mm:
                data += line
                continue
            mm = re.match(r"^send:", line)
            if mm:
                data += line
                continue
            mm = re.match(r"^reply:", line)
            if mm:
                data += line
                continue
            mm = re.match(r"^header:", line)
            if mm:
                data += line
                continue
            f.seek(last_tell)
            break
        return data
