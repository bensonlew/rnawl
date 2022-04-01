# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .admin import Admin
import sys
import os
import re
import time
import datetime
import json
import urlparse
import urllib
reload(sys)
sys.setdefaultencoding('utf8')


class ApilogAction(Admin):
    def __init__(self):
        super(ApilogAction, self).__init__()
        self.nav_index = 4
        self.title = "ApiLog日志"
        self.check_required(["workflow_id"])
        self.workflow_id = None
        self.start_time = None
        self.end_time = None
        self.log_files = []
        self.time = None
        self.area = None

    def GET(self):
        model = self.get_model("workflow")
        data = model.get_by_workflow_id(self.params.workflow_id)
        if data:
            self.workflow_id = data.workflow_id
            if data.has_run:
                self.log_files = self.get_log_files(data)
            else:
                return self.failed("尚未开始运行，无法查看日志！")
        else:
            return self.failed("workflow不存在，无法查看！", "workflow")
        return self.render.admin.apilog(self)

    # def POST(self):
    #     model = self.get_model("workflow")
    #     data = model.get_by_workflow_id(self.params.workflow_id)
    #     self.log_files = self.get_log_files(data)
    #     self.check_required(["time", "area"])
    #     self.time = time.mktime(time.strptime(self.params.time, '%Y-%m-%d %H:%M:%S'))
    #     self.area = int(self.params.area)
    #     return self.code

    def get_log_files(self, data):
        if self.params.addtime:
            self.start_time = datetime.datetime.strptime(self.params.addtime,
                                                         '%Y-%m-%d %H:%M:%S') + datetime.timedelta(minutes=-1)
            self.end_time = datetime.datetime.strptime(self.params.addtime,
                                                       '%Y-%m-%d %H:%M:%S') + datetime.timedelta(minutes=1)
        else:
            self.start_time = data.run_time + datetime.timedelta(minutes=-2)
            if data.end_time:
                self.end_time = data.end_time + datetime.timedelta(minutes=1)
            else:
                self.end_time = datetime.datetime.now()
        date_list = []
        for n in range((self.end_time - self.start_time).days+1):
            x = (self.start_time + datetime.timedelta(days=n)).strftime('%Y%m%d')
            date_list.append(x)
        log_dir = self.config.UPDATE_LOG
        run_host = re.sub("\.local$", "", data.server)
        return ["%s.%s.log" % (os.path.join(log_dir, x), run_host) for x in date_list]

    @property
    def code(self):
        in_area = False
        start_time = time.mktime(self.start_time.timetuple())
        end_time = time.mktime(self.end_time.timetuple())
        data = ""
        work_dir = os.environ['HOME']
        find_send = False
        trace_start = False
        for log in self.log_files:
            if os.path.exists(log):
                with open(log, "r") as f:
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        m = re.match(r'^(\d{4}\-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', line)
                        if m:
                            line_time = time.mktime(time.strptime(m.group(1), '%Y-%m-%d %H:%M:%S'))
                            # if self.time:
                            #     if abs(self.time - line_time) > self.area:
                            #         in_area = False
                            #         continue
                            #     else:
                            #         in_area = True
                            if (line_time >= start_time) and (line_time <= end_time):
                                in_area = True
                            else:
                                continue
                        if in_area:
                            if trace_start and re.match(r"^\w+:", line):
                                data += line
                                trace_start = False
                            if re.match(r"^Traceback", line):
                                trace_start = True
                            if trace_start:
                                m = re.match(r"^\s+File \"(.+)\", line (\d+), (.*)$", line)
                                if m:
                                    link_str = "<a target=\"_blank\" href=\"pymodule?path=%s&line=%s\">%s</a>" % \
                                               (os.path.relpath(m.group(1), work_dir), m.group(2), m.group(1))
                                    data += "  File \"%s\", line %s, %s\n" % (link_str, m.group(2), m.group(3))
                                else:
                                    data += line
                                continue
                            mm = re.search(r"(\{.*\})", line)
                            if mm:
                                try:
                                    json_data = json.loads(mm.group(1))
                                except Exception:
                                    data += line
                                    continue
                                else:
                                    if "id" in json_data.keys():
                                        if json_data["id"] != self.workflow_id:
                                            continue
                                    elif "task_id" in json_data.keys():
                                        if json_data["task_id"] != self.workflow_id:
                                            continue
                                        else:
                                            if "files" in json_data["data"]["sync_task_log"].keys():
                                                find_send = True
                                    else:
                                        if line.find("%s:" % self.workflow_id) == -1:
                                            continue
                                    if len(mm.group(1)) > 200:
                                        button = '<button type="button" class="btn btn-sm btn-default json" data-json="'\
                                                 + urllib.quote(mm.group(1)) + '" data-toggle="modal" ' \
                                                                             'data-target="#ApiLogJsonModal" ' \
                                                                             'aria-label="Json数据">JSON</button>'
                                        data += line.replace(mm.group(1), button)
                                    else:
                                        data += line
                                    continue
                            mmm = re.search(r"^send:.*&(sync_task_log=.*)\'", line)
                            if mmm:
                                json_data = json.loads(urlparse.parse_qs(mmm.group(1))["sync_task_log"][0])

                                if "task" in json_data.keys() and json_data["task"]["task_id"] != self.workflow_id:
                                    if find_send:
                                        if self.workflow_id.find(json_data["task"]["task_id"]) != 0:
                                            self.read_send(f)
                                            continue
                                        else:
                                            find_send = False
                                    else:
                                        self.read_send(f)
                                        continue
                                button = '<button type="button" class="btn btn-sm btn-default" data-urlcode="' \
                                         + mmm.group(1) + '" data-toggle="modal" ' \
                                         'data-target="#urlModal" aria-label="URLENCODE数据">URI CODE</button>'
                                data += line.replace(mmm.group(1), button)
                                data += self.read_send(f)
                                continue
                            if line.find("%s" % self.workflow_id) == -1:
                                continue
                            data += line
        return data

    @staticmethod
    def read_send(f):
        data = ""
        while True:
            line = f.readline()
            data += line
            mm = re.search(r"^Return page:", line)
            if mm:
                line = f.readline()
                data += line
                break
        return data
