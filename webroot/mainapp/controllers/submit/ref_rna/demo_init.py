# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from biocluster.wpm.client import worker_client, wait


class DemoInitAction(object):
    """
    demo设置
    """
    def __init__(self):
        super(DemoInitAction, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        requires = ['type', 'task_id', 'setup_type']
        for i in requires:
            if not (hasattr(data, i)):
                return json.dumps({"success": False, "info": "缺少%s参数!" % i})
        try:
            demo_number = data.demo_number
        except:
            demo_number = 2
        workflow_id = "DemoInit_" + "{}_{}".format(data.task_id, datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        if data.type == "ref_rna":
            data = {
              'id': workflow_id,
              'stat_id': 0,
              "type": "workflow",
              'name': "copy_demo.demo_init",  # 需要配置
              'client': data.client,
              "IMPORT_REPORT_DATA": False,
              "IMPORT_REPORT_AFTER_END": False,
              'options': {
                  "task_id": data.task_id,
                  "type": data.type,
                  "setup_type": data.setup_type,
                  "demo_number": demo_number
              }
            }
        try:
            worker = worker_client()
            info = worker.add_task(data)
            #print info
            if "success" in info.keys() and info["success"]:
                return {"success": True, "info": "任务提交成功%s" % (info["info"])}
            else:
                return {"success": False, "info": "任务提交失败%s" % (info["info"])}
        except:
            return {"success": False, "info": "任务提交失败"}
