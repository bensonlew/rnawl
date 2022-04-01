# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import web
import json
import datetime
from mainapp.models.mongo.submit.paternity_test_mongo import PaternityTest as PT
from mainapp.controllers.project.pt_controller import PtController
from mainapp.libs.param_pack import *
from mainapp.libs.signature import check_sig
from bson import ObjectId


class PtDatasplit(PtController):

    def __init__(self):
        super(PtDatasplit, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # print(data)
        # client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        # params_name = ['message_table', 'data_dir', 'family_table', 'submit_location', 'member_id']
        # params_name = ['message_table', 'data_dir', 'family_table', 'member_id', 'customer_table']
        params_name = ['message_table', 'data_dir', 'member_id']
        for param in params_name:
            if not hasattr(data, param):
                info = {'success': False, 'info': '缺少{}参数'.format(param)}
                return json.dumps(info)
        if not hasattr(data, "family_table"):
            if not hasattr(data, "customer_table"):
                info = {'success': False, 'info': '缺少客户信息表参数，无法进行分析'}
                return json.dumps(info)
        task_name = 'paternity_test.pt_datasplit'
        task_type = 'workflow'
        params_json = {
            'message_table': data.message_table,
            'data_dir': data.data_dir,
            # 'family_table': data.family_table,
            'member_id': data.member_id,
            # 'customer_table': data.customer_table,
        }
        if hasattr(data, "family_table"):
            params_json['family_table'] = data.family_table
        if hasattr(data, "customer_table"):
            params_json['customer_table'] = data.customer_table
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ('params', params),
            ('name', datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            ('desc', '亲子鉴定数据拆分'),
            ('member_id', data.member_id),
            ('status', 'start'),
            ('time', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_table_id = PT().insert_main_table('sg_pt_datasplit', mongo_data)
        update_info = {str(main_table_id): 'sg_pt_datasplit'}
        options = {
            "message_table": data.message_table,
            "data_dir": data.data_dir,
            # "family_table": data.family_table,
            "pt_data_split_id": str(main_table_id),
            "update_info": json.dumps(update_info),
            "member_id": data.member_id,
            # "customer_table": data.customer_table,
        }
        if hasattr(data, "family_table"):
            options['family_table'] = data.family_table
        if hasattr(data, "customer_table"):
            options['customer_table'] = data.customer_table
        sheet_data = self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params)
        # print "*********"
        # print sheet_data
        task_info = super(PtDatasplit, self).POST()
        return json.dumps(task_info)
