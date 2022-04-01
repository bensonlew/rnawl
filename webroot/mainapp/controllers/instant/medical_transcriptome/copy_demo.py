# -*- coding: utf-8 -*-
import web
import json
from mainapp.libs.signature import check_sig
from mainapp.controllers.core.basic import Basic
from biocluster.core.function import filter_error_info
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.core.base import Base
from biocluster.wpm.client import *
import random
import datetime
import unittest
import os
from mbio.packages.medical_transcriptome.copy_demo import CopyDemoMongo


class CopyDemoAction(Base):
    def __init__(self):
        super(CopyDemoAction, self).__init__()
        self.input_data = web.input()
        self._project_type = self.input_data.type

    @check_sig
    def POST(self):
        self.input_data = web.input()
        requires = [
            'type', 'task_id', 'target_task_id', 'target_project_sn',
            'target_member_id', 'target_member_type',
        ]
        for i in requires:
            if not (hasattr(self.input_data, i)):
                return json.dumps({"success": False, "info": "Lack argument: %s!" % i})
        old_task_id = self.input_data.task_id
        new_task_id = self.input_data.target_task_id
        new_project_sn = self.input_data.target_project_sn
        new_member_id = self.input_data.target_member_id
        project_type = self.input_data.type
        try:
            copy_demo = CopyDemoMongo(old_task_id, new_task_id, new_project_sn, new_member_id, project_type)
            copy_demo.run()
            run_info = dict(success=True)
            run_info['info'] = 'copy_demo_of_{}'.format(self.input_data.task_id)
            time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            col = self.db["sg_task"]
            demo_task_id = self.input_data.target_task_id
            run_info['demo_task_id'] = demo_task_id  # 返回给前端的run_info信息增加demo task id字段
            col.update({"task_id": demo_task_id}, {"$set": {'status': "end"}}, upsert=True)
            col.update({"task_id": demo_task_id}, {"$set": {"is_demo": 1}}, upsert=True)
            col.update({"task_id": demo_task_id}, {"$set": {"created_ts": time}})
            col.update({"task_id": demo_task_id}, {"$set": {"member_type": self.input_data.target_member_type}})
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
        cmd += "i/project_demo/copy_demo "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id='RefrnaV2_7320',
            target_task_id="copy_demo_ffff",
            type="ref_rna_v2",
            target_project_sn="test",
            target_member_id='test',
            target_member_type='1',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
