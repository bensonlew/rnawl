# -*- coding: utf-8 -*-
# __author__ = 'hongyu'
import web
import json
import datetime
from bson import SON
from mainapp.libs.param_pack import *
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.submit.med_mongo import MedMongo as PT
from mainapp.controllers.project.pt_controller import PtController


class PtDedupAction(PtController):
    """
    该接口用于全库排查守护进程
    """
    def __init__(self):
        super(PtDedupAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        for params_name in ["back_time", "err_min", "mem"]:
            if not (hasattr(data, params_name)):
                info = {"success": False, "info": '缺少参数' + params_name}
                return json.dumps(info)
        task_name = 'medical.paternity_test_v2.report.pt_daily_dedup'
        task_type = 'workflow'
        task_id = "Pt_Daily_Dedup_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        params_json = {
            "back_time": data.back_time,
            "err_min": data.err_min,
            "mem": data.mem,
            "time": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ('params', params),
            ('status', 'start'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_table_id = PT('pt_v2').add_main_id('sg_daily_dedup')
        update_info = {str(main_table_id): 'sg_daily_dedup'}
        options = {
            "back_time": data.back_time,
            "err_min": data.err_min,
            "mem": data.mem,
            "update_info": json.dumps(update_info),
            'main_table_data': SON(mongo_data)
        }
        self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params, db_type='pt_v2')
        task_info = super(PtDedupAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': task_id}}
        return json.dumps(task_info)
