# -*- coding: utf-8 -*-

import web
import json
import os
import unittest
from mainapp.controllers.project.delete_demo_controller import DeleteDemoController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.core.base import Base
from mainapp.controllers.core.basic import Basic
from biocluster.core.function import filter_error_info
from biocluster.wpm.client import *
import datetime

class DeleteDemoV2Action(DeleteDemoController):
    def __init__(self):
        super(DeleteDemoV2Action, self).__init__()
        self.input_data = web.input()
        self._project_type = self.input_data.project_type

    @check_sig
    def POST(self):
        print(self.input_data)
        requires = ['project_type', 'task_id', 'main_id']
        for i in requires:
            if not (hasattr(self.input_data, i)):
                return json.dumps({"success": False, "info": "Lack argument: %s!" % i})
        info = self.delete_demo.update_main_table('sg_task_delete', self.input_data.task_id, self.input_data.main_id, 'deleting')
        if not info['success']:
            return json.dumps({"success": False, "info": info['info']})
        workflow_id = "Delete_" + "{}_{}".format(self.input_data.task_id, datetime.datetime.now().strftime("%H%M%S%f")[:-3])
        data_json = {
            'id': workflow_id,
            'stage_id': 0,
            'name': "project_demo.delete_demo_v2",  # 需要配置
            'type': 'workflow',  # 可以配置
            'client': self.input_data.client,
            "IMPORT_REPORT_DATA": False,
            "IMPORT_REPORT_AFTER_END": False,
            'options': {
                "task_id": self.input_data.task_id,
                "project_type": self.input_data.project_type,
                "db_version": info['db_version'],
                "main_id": self.input_data.main_id
            }
        }

        workflow_client = Basic(data=data_json, instant=False)
        try:
            run_info = workflow_client.run()
            run_info['info'] = filter_error_info(run_info['info'])
            return json.dumps(run_info)
        except Exception as e:
            return json.dumps({"success": False, "info": "running error: %s" % filter_error_info(str(e))})

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/project_demo/delete_demo_v2 "
        cmd += "-b http://bcl.tsg.com "
        cmd += "-dbversion 0 "
        args = dict(
            task_id="8cre_s5076tktj23b3719lppvqu",
            project_type="medical_transcriptome",
            main_id="60c2d6eda6dde5652a06a694"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
