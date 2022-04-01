# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

import web
import json
import datetime
from bson import ObjectId
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.arghub_controller import ArghubController


class ArghubAction(ArghubController):
    def __init__(self):
        super(ArghubAction, self).__init__(instant=False)

    def POST(self):
        data = web.input()
        print data
        if not hasattr(data, "input_type"):
            info = {"success": False, "info": "PARAMETERS MISSING: input_type" }
            return json.dumps(info)
        input_type = data.input_type
        if input_type == "read":
            if not(hasattr(data, "read1") and hasattr(data, "read2")):
                info = {"success": False, "info": "input_type 为 read是必须配置 read1 read2参数"}
                return json.dumps(info)
        else:
            if not hasattr(data, input_type) and not hasattr(data, "sequence"):
                info = {"success": False, "info": "input_type 为 {0} 是必须配置 {0}参数".format(input_type)}
                return json.dumps(info)

        task_name = 'arghub.arghub'
        module_type = 'workflow'
        task_type = 2
        params_json = {'input_type': data.input_type}
        if data.input_type == "read":
            params_json["read1"] = data.read1
            params_json["read2"] = data.read2
        else:
            if hasattr(data, input_type):
                params_json[data.input_type] = data[data.input_type]
            else:
                params_json["sequence"] = data["sequence"]
        for a in ["gff", "aligner", "evalue", "identity"]:
            if hasattr(data, a):
                params_json[a] = data[a]
        mongo_data = [
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', 'processing'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('input_type', data.input_type)
        ]
        options = params_json.copy()
        main_table_name = "Arghub_{}".format(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.arghub.insert_main_table('analysis', mongo_data)
        update_info = {str(main_table_id): 'analysis'}
        options['main_id'] = str(main_table_id)
        options['update_info'] = json.dumps(update_info)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=data['task_id'], mem_id=data.member_id,
                            module_type=module_type, params=params_json)
        task_info = super(ArghubAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

