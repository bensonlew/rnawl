# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last modify 20180103
import web
import json
import datetime
from mainapp.models.mongo.submit.med_mongo import MedMongo as Nipt
from mainapp.controllers.project.nipt_controller import NiptController
from bson import ObjectId
from mainapp.libs.signature import check_sig
from bson import SON


class NiptAnalysisV2Action(NiptController):
    def __init__(self):
        super(NiptAnalysisV2Action, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['main_id', 'bw', 'bs', 'ref_group']
        for param in params_name:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!!" % param}
                return json.dumps(info)
        task_info = Nipt("nipt_v2").get_one("sg_main", "_id", ObjectId(data.main_id))
        if not task_info:
            info = {"success": False, "info": "不存在 %s 的主表，请确认参数是否正确！!"% data.main_id}
            return json.dumps(info)
        task_name = 'medical.nipt_v2.report.nipt_analysis'
        task_type = 'workflow'
        main_table_name = 'bw-' + data.bw + '_bs-' + data.bs + '_ref-' + data.ref_group + \
                          datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")[:-3]
        params_json = {
            'bw': data.bw,
            'bs': data.bs,
            'ref_group': int(data.ref_group),
            'main_id': data.main_id
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ('batch_id', task_info['batch_id']),
            ('sample_id', task_info['sample_id']),
            ('status', 'start'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('params', params),
            ('name', main_table_name),
            ("main_id", ObjectId(data.main_id))
        ]
        main_table_id = Nipt('nipt_v2').add_main_id('sg_interaction')
        update_info = {str(main_table_id): 'sg_interaction'}
        options = {
            "bed_file": task_info['sample_id'],
            "bw": data.bw,
            "bs": data.bs,
            "update_info": json.dumps(update_info),
            "ref_group": data.ref_group,
            "nipt_task_id": str(main_table_id),
            'main_table_data': SON(mongo_data)
        }
        self.set_sheet_data(name=task_name, options=options, module_type=task_type, params=params, db_type='nipt_v2')
        task_info = super(NiptAnalysisV2Action, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
