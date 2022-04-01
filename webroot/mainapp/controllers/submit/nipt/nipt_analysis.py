# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last modify 20170522
import web
import json
import datetime
from mainapp.models.mongo.submit.nipt_mongo import NiptMongo as Nipt
from mainapp.controllers.project.nipt_controller import NiptController
from bson import ObjectId
from mainapp.libs.signature import check_sig
from bson import SON


class NiptAnalysis(NiptController):
    def __init__(self):
        super(NiptAnalysis, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['main_id', 'bw', 'bs', 'ref_group']
        for param in params_name:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!!" % param}
                return json.dumps(info)
        task_info = Nipt().get_sample_info(data.main_id)
        if not task_info:
            info = {"success": False, "info": "不存在 %s 的主表，请确认参数是否正确！!"%(data.main_id)}
            return json.dumps(info)
        task_name = 'nipt.report.nipt_analysis'
        task_type = 'workflow'
        # main_table_name = 'Nipt_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        main_table_name = 'bw-' + data.bw + '_bs-' + data.bs + '_ref-' + data.ref_group
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
        # main_table_id = Nipt().insert_none_table('sg_interaction')
        main_table_id = Nipt().insert_main_table('sg_interaction', mongo_data)
        update_info = {str(main_table_id): 'sg_interaction'}
        options = {
            "bed_file": task_info['sample_id'],
            "bw": data.bw,
            "bs": data.bs,
            "update_info": json.dumps(update_info),
            "ref_group": data.ref_group,
            "nipt_task_id": str(main_table_id),
            # 'main_table_data': SON(mongo_data)
        }
        self.set_sheet_data(name=task_name, options=options, module_type=task_type, params=params)
        task_info = super(NiptAnalysis, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
