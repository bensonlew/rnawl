# -*- coding: utf-8 -*-
# __author__ = 'xuting'

import web
import json
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow


class Datasplit(object):
    """
    数据拆分接受接口
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        if not hasattr(data, "content") or not hasattr(data, "task_id"):
            info = {"success": False, "info": "Json内容不正确!!"}
            print json.dumps(info)
            return json.dumps(info)
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        json_obj = dict()
        json_obj["options"] = dict()
        json_obj["options"]["sample_info"] = json.loads(data.content)
        json_obj["name"] = "datasplit.datasplit"
        json_obj["USE_DB"] = True
        json_obj["UPDATE_STATUS_API"] = "split_data"
        json_obj['client'] = client
        json_obj['type'] = "workflow"
        json_obj['id'] = data.task_id
        json_obj["stage_id"] = data.task_id
        json_obj["to_file"] = ["datasplit.export_sample_info(sample_info)"]
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(json_obj['id'])
        if len(workflow_data) > 0:
            info = {"success": False, "info": "pipeline流程ID重复!"}
            print json.dumps(info)
            return json.dumps(info)
        else:
            insert_data = {"client": client,
                           "workflow_id": json_obj['id'],
                           "json": json.dumps(json_obj),
                           "ip": web.ctx.ip}
            workflow_module.add_record(insert_data)
            info = {"success": True, "info": "添加队列成功!"}
            print json.dumps(info)
            return json.dumps(info)
