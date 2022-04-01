# -*- coding: utf-8 -*-
# __author__ = "zouguanqing"
import web
import json
import datetime
from mainapp.controllers.project.metagenomic_controller import MetagenomicController

from mainapp.libs.signature import check_sig


class FunctionSetAction(MetagenomicController):
    def __init__(self):
        super(FunctionSetAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['is_ref', "anno_type", "level", "level_list", 'task_type', 'task_id','name', 'submit_location', "type"] #type: env or medical

        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" ,"variables": [argu], "code" : "C2403701"}
                return json.dumps(info)
        task_name = "metagenomic.report.function_set"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type
        task_id = data.task_id
        geneset_table = self.metagenomic.find_origin_geneset_info(data.task_id)
        project_sn = geneset_table["project_sn"]

        name =data.name
        params_json = {
            "database": data.anno_type,
            "submit_location": data.submit_location,
            "task_type": int(task_type),
            "type" : data.type,

            #"level": data.level,
            #"name":name
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("name", name),
            ("desc", ""),
            ("func_desc",""),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("anno_type" , data.anno_type),
            ("level" , data.level),
            ("member", data.level_list),
            ("submit_location" ,data.submit_location),
            ("is_ref", int(data.is_ref)),
            ("type",data.type)
        ]

        options = {
            "database" : data.anno_type,
            "level" : data.level,
            "member" : data.level_list,
        }

        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))

        #self.from_id_get_result('func_set',data.func_id)
        if hasattr(data, 'funcset_id') and  data.funcset_id != '':
            options["main_id"] = str(data.funcset_id)
        else:
            main_table_id = self.metagenomic.insert_main_table("func_set", mongo_data)
            options["main_id"] = str(main_table_id)

        update_info = {options["main_id"] : "func_set"}
        options["update_info"] = json.dumps(update_info)

        self.set_sheet_data(name=task_name, options=options, main_table_name='Function_Set',
                            module_type=module_type, project_sn=project_sn,
                            task_id=task_id, params=params_json)
        task_info = super(FunctionSetAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': options["main_id"], 'name': 'Function_Set'}}
        return json.dumps(task_info)

