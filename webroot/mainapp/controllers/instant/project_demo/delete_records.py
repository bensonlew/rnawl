# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

import datetime
import json
import os
import unittest

import web
from biocluster.core.function import filter_error_info
from biocluster.wpm.client import *

from mainapp.controllers.core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.core.base import Base


class DeleteRecordsAction(Base):
    def __init__(self):
        super(DeleteRecordsAction, self).__init__()
        self.input_data = web.input()
        self._project_type = self.input_data.type

    @check_sig
    def POST(self):
        requires = ['type', 'task_id', 'submit_location', 'status']
        for i in requires:
            if not (hasattr(self.input_data, i)):
                return json.dumps({"success": False, "info": "Lack argument: %s!" % i})

        workflow_id = "Delete_" + "{}_{}".format(self.input_data.task_id,
                                                 datetime.datetime.now().strftime("%H%M%S%f")[:-3])
        data_json = {
            'id': workflow_id,
            'stage_id': 0,
            'name': "project_demo.delete_records",  # 需要配置
            'type': 'workflow',  # 可以配置
            'client': self.input_data.client,
            # 'project_sn': data.project_sn,
            "IMPORT_REPORT_DATA": False,
            "IMPORT_REPORT_AFTER_END": False,
            'options': {
                "task_id": self.input_data.task_id,
                "project_type": self.input_data.type,
                "submit_location": self.input_data.submit_location,
                "status": self.input_data.status
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
        cmd = 'python /mnt/ilustre/users/isanger/sanger_bioinfo/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/project_demo/delete_records "
        cmd += "-b http://bcl.i-sanger.com "
        args = dict(
            task_id="majorbio_283598",
            type="denovo_rna_v2",
            submit_location="diff_detail",
            status="failed"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
