
# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.models.mongo.meta import Meta
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig



class VpaAction(MetaController):
    def __init__(self):
        super(VpaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = [ "otu_id",'submit_location', 'group_detail','task_type',
                         'group_id', 'env_id','env_group_id', 'env_detail', 'level_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables":[ argu], "code" : "C2204101"}
                return json.dumps(info)

        task_name = 'meta.report.vpa'
        module_type = 'workflow'
        # task_type = 2
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", "code" : "C2204102"}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'level_id': data.level_id,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'env_id': data.env_id,
            'env_detail': group_detail_sort(data.env_detail),
            'env_group_id': data.env_group_id,
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            #'env_labs' : data.env_labs
        }
        data.env_detail = self.meta.get_new_env_detail(data.env_detail, data.env_id)
        print(data)

        main_table_name = 'VPA_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.meta.insert_main_table('sg_vpa', mongo_data)
        update_info = {str(main_table_id): 'sg_vpa'}
        options = {
            'otu_table': data.otu_id,
            'otu_id': data.otu_id,
            'level': int(data.level_id),
            'env_id': data.env_id,
            'env_file': data.env_id,
            'env_labs': ','.join([','.join(i) for i in json.loads(data.env_detail).values()]),
            'env_detail': data.env_detail,
            'env_group':data.env_detail,
            # 'group_table': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'group_id': data.group_id,
            'group': data.group_id
        }

        to_file = ["meta.export_otu_table_by_level(otu_table)",
                    "meta.export_group_table_by_detail(group)", "env.export_float_env_regression(env_file)"]
        to_file.append('meta.export_env_table_by_detail(env_group)')

        self.set_sheet_data(name=task_name, options=options, main_table_name="VPA/"+main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(VpaAction, self).POST()
        task_info['content'] = {
            'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)