# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last modify 20170821
import web
import json
import datetime
from mainapp.models.mongo.submit.nipt_mongo import NiptMongo as Nipt
from mainapp.controllers.project.nipt_controller import NiptController
from bson import ObjectId
from bson import SON
from mainapp.libs.signature import check_sig
import os


class MakeRef(NiptController):
    def __init__(self):
        super(MakeRef, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # print data
        params_name = ['sample_txt', 'ref_group']
        for param in params_name:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!!" % param}
                return json.dumps(info)
        task_name = 'nipt.report.make_ref'
        task_type = 'workflow'
        if str(data.client) == "client03":
            sample_file = os.path.join('/mnt/ilustre/tsanger-data', data.sample_txt)
        elif str(data.client) == "client01":
            sample_file = os.path.join('/mnt/ilustre/data/', data.sample_txt)
        else:
            return json.dumps({"success": False, 'info': "client不合法，client必须是client01 or client03"})
        print sample_file
        main_table_name = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3] + '_ref-' + data.ref_group
        params_json = {
            'sample_txt': data.sample_txt,
            'ref_group': int(data.ref_group),
            'file_id': str(data.file_id),
            'sample_name': data.sample_name
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ('sample_txt', data.sample_name),
            ('status', 'start'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('params', params),
            ('name', main_table_name),
            ("use_status", "当前在使用")
        ]
        main_table_id = Nipt().insert_main_table('sg_ref_main', mongo_data)
        update_info = {str(main_table_id): 'sg_ref_main'}

        options = {
            "sample_txt": sample_file,
            "update_info": json.dumps(update_info),
            "ref_group": data.ref_group,
            "nipt_task_id": str(main_table_id)
        }
        self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params)
        task_info = super(MakeRef, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
