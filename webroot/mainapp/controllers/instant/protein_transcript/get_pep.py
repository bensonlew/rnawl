# -*- coding: utf-8 -*-
import web
import json
import os
import unittest
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.core.base import Base
from mainapp.controllers.core.basic import Basic
from biocluster.core.function import filter_error_info
from biocluster.wpm.client import *
import datetime
# from mbio.packages.project_demo.delete_demo import DeleteRelaMongo
# import sys


class GetPepAction(Base):
    def __init__(self):
        super(GetPepAction, self).__init__()
        self.input_data = web.input()
        self._project_type = self.input_data.type

    @check_sig
    def POST(self):
        requires = ['type', 'task_id', 'dir_suffix']
        for i in requires:
            if not (hasattr(self.input_data, i)):
                return json.dumps({"success": False, "info": "Lack argument: %s!" % i})

        workflow_id = "{}_{}".format(self.input_data.task_id, datetime.datetime.now().strftime("%H%M%S%f")[:-3])
        data_json = {
            'id': workflow_id,
            'stage_id': 0,
            'name': "protein_transcript.report.get_pep",  # 需要配置
            'type': 'workflow',  # 可以配置
            'client': self.input_data.client,
            # 'project_sn': data.project_sn,
            "IMPORT_REPORT_DATA": False,
            "IMPORT_REPORT_AFTER_END": False,
            'options': {
                "task_id": self.input_data.task_id,
                "project_type": self.input_data.type,
                "dir_suffix": self.input_data.dir_suffix
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
        cmd = 'python /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/protein_transcript/get_pep "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            task_id="onr1_m8cnapnpg3c8bu276m2o3c",
            task_type="1",
            dir_suffix="12345",
            type="ref_rna_v2",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
