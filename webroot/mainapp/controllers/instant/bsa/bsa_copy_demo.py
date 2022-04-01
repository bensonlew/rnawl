# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.0305

import web
import json
import random
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.controllers.core.basic import Basic
# # from biocluster.wpm.client import worker_client, wait
from biocluster.core.function import filter_error_info
from mainapp.controllers.project.bsa_controller import BsaController
import datetime


class BsaCopyDemoAction(BsaController):
    """
    BSA demo拉取接口
    """
    def __init__(self):
        super(BsaCopyDemoAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        # params = ["type", "task_id", "target_task_id", "target_project_sn", "target_member_id", "target_member_type", "target_cmd_id"]
        params = ["type", "task_id", "target_task_id", "target_project_sn", "target_member_id"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" , "variables": [param], "code" : "C1200105"}
                return json.dumps(info)
        workflow_id = self.get_new_id(data.task_id)
        if data.type == "bsa":
            project_type = "dna_bsa"
        elif data.type == "dna_gmap.gmap":
            project_type = "dna_gmap"
        elif data.type == "dna_evolution.evolution":
            project_type = "dna_evolution"
        elif data.type == "noref_wgs.noref_wgs":
            project_type = "noref_wgs"
        elif data.type == "wgs_v2.wgs_v2":
            project_type = "dna_wgs_v2"
        else:
            info = {"success": False, "info": "项目类型: %s不对，请检查!" , "variables": [data.type], "code" : "C1200106"}
            return json.dumps(info)
        data_json = {
            "id": workflow_id,
            "stage_id": 0,        
            "name": "copy_demo.bsa_copy_demo",
            "type": "workflow",
            "client": data.client,
            "project_sn": data.target_project_sn,
            "options": {
                "task_id": data.task_id,
                "target_task_id": data.target_task_id,
                "target_project_sn": data.target_project_sn,
                "target_member_id": data.target_member_id,
                "project_type": project_type
            }
        }
        worker_client = Basic(data=data_json, instant=True)
        try:
            run_info = worker_client.run()
            run_info["info"] = filter_error_info(run_info["info"])
            return json.dumps(run_info)
        except Exception as e:
            var = []
            var.append(filter_error_info[str(e)])
            # return json.dumps({"success": False, "info": "拉取demo运行出错:{}".format(filter_error_info[str(e)])})
            return json.dumps({"success": False, "info": "拉取demo运行出错:%s" % filter_error_info[str(e)],
                               "code": "C1200103", "variables": var})

    def get_new_id(self, task_id):
        """
        根据旧的ID生成新的workflowID，固定为旧的后面用“_”，添加两次随机数
        """
        new_id = "{}_{}_{}".format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id)
        return new_id
