# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0510

import web
import json
import random
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.controllers.core.basic import Basic
# # from biocluster.wpm.client import worker_client, wait
from biocluster.core.function import filter_error_info
# import class DnaController
from mainapp.controllers.project.dna_controller import DnaController
import datetime


class WgsCopyDemoAction(DnaController):
    """
    WGS demo拉取接口
    """

    def __init__(self):
        super(WgsCopyDemoAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data_json = {}
        data = web.input()
        print "=================="
        print data
        params = [
            "type",
            "task_id",
            "target_task_id",
            "target_project_sn",
            "target_member_id",
            "target_member_type",
            "target_cmd_id"]         # target_member_type项目类型的权限，是生信还是其他
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {
                    "success": False,
                    "info": "缺少%s参数!" %
                    param,
                    "code": "C3101601",
                    "variables": var}
                return json.dumps(info)
        workflow_id = self.get_new_id(data.task_id)
        if data.type == "wgs.wgs":
            data_json = {
                "id": workflow_id,
                "stage_id": 0,
                "name": "copy_demo.wgs_copy_demo",
                "type": "workflow",
                "client": data.client,
                "project_sn": data.target_project_sn,
                "options": {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id
                }
            }
        else:
            info = {
                "success": False,
                "info": "项目类型不对，请检查",
                "code": "C3101602",
                "variables": ""}
            return json.dumps(info)
        worker_client = Basic(data=data_json, instant=True)
        try:
            run_info = worker_client.run()
            run_info["info"] = filter_error_info(run_info["info"])
            return json.dumps(run_info)
        except Exception as e:
            var = []
            var.append(filter_error_info[str(e)])
            return json.dumps({"success": False, "info": "拉取demo运行出错:%s" % (
                filter_error_info[str(e)]), "code": "C3101603", "variables": var})

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
