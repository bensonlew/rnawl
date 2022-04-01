# -*- coding: utf-8 -*-
# __author__ = 'HD'
# create 20210429

import web
import json
from mainapp.libs.signature import check_sig
from biocluster.wpm.client import worker_client


class S3filedeleteAction(object):

    @check_sig
    def POST(self):
        data = web.input()
        for i in ["project_sn", "task_id", "unique_id"]:
            if not hasattr(data, i):
                msg = {"success": False, "info": "缺少参数:%s" % i}
                return json.dumps(msg)
            if getattr(data, i).strip() == "" and i != "unique_id":
                msg = {"success": False, "info": "参数%s不能为空" % i}
                return json.dumps(msg)
        json_obj = {
            "project_id": data.project_sn,
            "task_id": data.task_id,
            "unique_id": data.unique_id
        }
        if hasattr(data, "wfm_port"):
            json_obj["wfm_port"] = data.wfm_port
        response = worker_client().add_filedelete_task(json_obj)
        if response:
            if "success" in response.keys() and response["success"]:
                info = {"success": True, "info": "添加任务到队列成功!", "workflow_id": json_obj['task_id']}
                return json.dumps(info)
        else:
            if response:
                resinfo = response["info"]
            else:
                resinfo = "请求GRPC失败"
            info = {"success": False, "info": "添加任务到队列失败：{}!".format(resinfo),
                    "workflow_id": json_obj['task_id']}
            return json.dumps(info)
