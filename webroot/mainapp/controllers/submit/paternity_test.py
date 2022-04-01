# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'

import web
import json
import datetime
import random
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.paternity_test_mongo import PaternityTest as import_mongodb
from mainapp.libs.param_pack import *

class PaternityTest(object):
    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        params_name = ['err_min', 'dedup_num','task']
        success = []
        if isinstance(int(data.dedup_num), int) == False:
            info = {"success": False, "info": "查重样本数量必须为整数"}
            return json.dumps(info)

        task_info = import_mongodb().get_query_info(data.task)
        flow_id = self.get_new_id(data.task)
        member_id = 'zml'
        (output_dir, update_api) = GetUploadInfo_pt(client, member_id, task_info['project_sn'], data.task,
                                                 'pt_report')
        json_data = {
	        "id": str(flow_id),
	        "project_sn": task_info['project_sn'],
	        # "task_id": data.task,
	        "type": "workflow",
	        "name": "paternity_test.report.pt_report",
	        "client": client,
	        "USE_DB": True,
	        "IMPORT_REPORT_DATA": True,
	        "UPDATE_STATUS_API": update_api,
	        "IMPORT_REPORT_AFTER_END": True,
	        "output": output_dir,
	        "options": {
		        "ref_fasta": str(task_info['ref_fasta']),
		        "targets_bedfile": str(task_info['targets_bedfile']),

		        "dad_id": str(task_info['dad_id']),
		        "mom_id": task_info['mom_id'],
		        "preg_id": task_info['preg_id'],
		        "ref_point": str(task_info['ref_point']),

		        "err_min": int(data.err_min),
		        "dedup_num": int(data.dedup_num),
		        "flow_id":flow_id
	        }
        }
        insert_data = {"client": client,
                           "workflow_id": flow_id,
                           "json": json.dumps(json_data),
                           "ip": web.ctx.ip
                           }
        workflow_module = Workflow()
        workflow_module.add_record(insert_data)
        info = {"success": True, "info": "提交成功!"}
        return json.dumps(info)

    def get_new_id(self, task_id):
        new_id = "%s_%s_%s" % (task_id, random.randint(1, 10000), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id)
        return new_id


# python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post paternity_test -c client03 -b http://192.168.12.102:8091 -n "err_min;dedup_num;task" -d "3;50;zml01"