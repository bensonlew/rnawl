# -*- coding: utf-8 -*-
import web
import json
import os
import unittest
import datetime
from mainapp.libs.signature import check_sig
from mainapp.controllers.core.basic import Basic
from biocluster.core.function import filter_error_info
from mainapp.models.mongo.core.base import Base
from biocluster.wpm.client import *


class DeleteDemoAction(Base):
    def __init__(self):
        super(DeleteDemoAction, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        requires = ['type', 'task_id']
        for i in requires:
            if not (hasattr(data, i)):
                return json.dumps({"success": False, "info": "缺少%s参数!" % i})
        workflow_id = "Delete_" + "{}_{}".format(data.task_id, datetime.datetime.now().strftime("%H%M%S%f")[:-3])

        if data.type == "dia_v3":
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "dia_v3.delete_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                # 'project_sn': data.project_sn,
                "IMPORT_REPORT_DATA": False,
                "IMPORT_REPORT_AFTER_END": False,
                'options': {
                    "task_id": data.task_id,
                    "project_type": data.type,
                }
            }
        else:
            info = {"success": False, "info": "项目类型type：{}不合法! 请输入合法的type,比如:dia_v3".format(data.type)}
            return json.dumps(info)

        workflow_client = Basic(data=data_json, instant=True)
        try:
            run_info = workflow_client.run()
            run_info['info'] = filter_error_info(run_info['info'])
            return json.dumps(run_info)
        except Exception as e:
            return json.dumps({"success": False, "info": "运行出错: %s" % filter_error_info(str(e))})


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/labelfree/delete_demo "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="copy_demo_test",
            type="labelfree",
            project_sn="labelfree",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
