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


class CopyDemoAction(Base):
    def __init__(self):
        super(CopyDemoAction, self).__init__()
        self._project_type = 'ref_rna_v2'

    @check_sig
    def POST(self):
        data = web.input()
        requires = [
            'type', 'task_id', 'target_task_id', 'target_project_sn',
            'target_member_id', 'target_member_type'
        ]
        for i in requires:
            if not (hasattr(data, i)):
                return json.dumps({"success": False, "info": "缺少%s参数!" % i})
        workflow_id = self.get_new_id(data.task_id)

        if data.type == "ref_rna_v2":
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "ref_rna_v2.copy_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "project_type": data.type,
                }
            }
        else:
            info = {"success": False, "info": "项目类型type：{}不合法!".format(data.type)}
            return json.dumps(info)

        workflow_client = Basic(data=data_json, instant=True)
        try:
            run_info = workflow_client.run()
            run_info['info'] = filter_error_info(run_info['info'])
            if data.type == "ref_rna_v2":
                time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                col = self.db["sg_task"]
                demo_task_id = data.target_task_id
                run_info['demo_task_id'] = demo_task_id  # 返回给前端的run_info信息增加demo task id字段
                col.update({"task_id": demo_task_id},{"$set": {'status': "end"}}, upsert=True)
                col.update({"task_id": demo_task_id}, {"$set": {"is_demo": 1}}, upsert=True)
                col.update({"task_id": demo_task_id}, {"$set": {"created_ts": time}})
                col.update({"task_id": demo_task_id}, {"$set": {"member_type": data.target_member_type}})
            return json.dumps(run_info)
        except Exception as e:
            return json.dumps({"success": False, "info": "运行出错: %s" % filter_error_info(str(e))})

    def get_new_id(self, task_id, otu_id=None):
        """
        根据旧的ID生成新的workflowID，固定为旧的后面用“_”，添加两次随机数或者一次otu_id一次随机数
        """
        if otu_id:
            new_id = "{}_{}_{}".format(task_id, otu_id[-4:], random.randint(1, 10000))
        else:
            new_id = "{}_{}_{}".format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"),
                                       random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id, otu_id)
        return new_id


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/ref_rna_v2/copy_demo "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id='tsg_28999',
            target_task_id="copy_demo_test",
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
