# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from biocluster.wpm.client import worker_client, wait


PACKAGE_URL = "demo_init"
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
        print data
        for i in requires:
            if not (hasattr(data, i)):
                return json.dumps({"success": False, "info": "缺少%s参数!" % i})
        try:
            demo_number = data.demo_number
        except:
            demo_number = 30
        if data.task_id == "" or data.task_id == " ":
            return json.dumps({"success": False, "info": "参数task_id不能为空!"})
        workflow_id = "DemoInit_" + "{}_{}".format(data.task_id, datetime.datetime.now().strftime("%H%M%S%f")[:-3])
        # if data.type == "ref_rna":
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
            print info
            if "success" in info.keys() and info["success"]:
                return {"success": True, "info": "demo备份任务提交成功,请两个小时后进行此demo的拉取或取消demo设置操作"}
            else:
                return {"success": False, "info": "demo备份任务提交失败,请重新设置"}
        except:
            return {"success": False, "info": "demo备份任务提交失败,请重新设置"}
