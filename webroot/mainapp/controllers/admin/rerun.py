# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .admin import Admin
import json
import re
from biocluster.wpm.client import worker_client
import os


class RerunAction(Admin):

    def __init__(self):
        super(RerunAction, self).__init__()
        self._table = "workflow"

    def POST(self):
        self.check_required(["wid", "step_names"])
        if self.params.step_names == "" and self.params.skip_all != "1":
            return json.dumps({"success": False, "info": "必须填写步骤名或指定跳过的Tool！"})
        model = self.get_model("workflow")
        data = model.get_by_id(self.params.wid)
        if len(data) > 0:
            json_data = json.loads(data.json)
            json_data["rerun"] = True
            json_data["run_time"] = data.run_time

            if data.work_dir:
                json_data["work_dir"] = data.work_dir
            if self.params.skip_all == "1":
                # tool_ids = self.get_ids_from_log(data.workflow_id, data.run_time, data.server)
                # tool_ids_2 = self.get_ids_from_db(data.workflow_id)
                # for i in tool_ids_2:
                #     if i not in tool_ids:
                #         tool_ids.append(i)
                # json_data["skip_tools"] = tool_ids
                json_data["skip_all_success"] = True
            if self.params.step_names:
                steps = re.split(r"\s*,\s*", self.params.step_names)
                json_data["skip_steps"] = steps
            if self.params.api_path == "None":
                json_data["UPDATE_STATUS_API"] = False
            elif self.params.api_path != "":
                json_data["UPDATE_STATUS_API"] = self.params.api_path
            if self.params.import_data == "1":
                json_data["IMPORT_REPORT_DATA"] = True
            else:
                json_data["IMPORT_REPORT_DATA"] = False
            if self.params.import_end == "1":
                json_data["IMPORT_REPORT_AFTER_END"] = True
            else:
                json_data["IMPORT_REPORT_AFTER_END"] = False
            try:
                worker = worker_client()
                info = worker.add_task(json_data)
                if "success" in info.keys() and info["success"]:
                    self.log("重运行workflow: %s " % data.workflow_id )
                    return json.dumps({"success": True, "info": "重运行成功！"})
                else:
                    return json.dumps({"success": False, "info": "任务投递失败: %s" % info["info"]})
            except Exception, e:
                return json.dumps({"success": False, "info": "WPM连接错误:%s" % e})

        else:
            return json.dumps({"success": False, "info": "没有找到对应的workflow！"})

    def get_ids_from_db(self, workflow_id):
        tool_model = self.get_model("tool")
        tools = tool_model.get_workflow_tools(workflow_id)
        tool_ids = []
        if tools:
            for t in tools:
                if t.run_id not in tool_ids:
                    tool_ids.append(t.run_id)
        return tool_ids

    def get_ids_from_log(self, workflow_id, run_time, server):
        timestr = run_time.strftime('%Y%m%d')
        log_dir = os.path.join(self.config.wpm_log_file, timestr)
        log_file = os.path.join(log_dir, "%s_%s.log" % (workflow_id, re.sub(r"\.local$", "", server)))
        tool_ids = []
        if os.path.exists(log_file):
            with open(log_file, "r") as f:
                for line in f:
                    line = line.strip('\n')
                    if not re.search(r"RPC", line):
                        continue
                    if not re.search(r"\'state\': \'finish\'", line):
                        continue
                    print line
                    m = re.search(r"\'id\': \'(\S+)\'", line)
                    if m:
                        wid = m.group(1)
                        if wid not in tool_ids:
                            tool_ids.append(wid)
        return tool_ids
