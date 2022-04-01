# -*- coding: utf-8 -*-
# __author__ = 'ysh'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class HomofilterAction(BacgenomeController):
    def __init__(self):
        super(HomofilterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:",data
        default_argu = ['task_id', 'all_samples','homology_id', 'filter_type','submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if data.filter_type not in ["core","unique","defined"]:
            info = {'success': False, 'info': 'filter must be core, filter_type or defined!'}
            return json.dumps(info)
        if not hasattr(data, "filter") and not hasattr(data, "select"):
            info = {'success': False, 'info': 'argu filter or select must be at last exists one!'}
            return json.dumps(info)
        task_name = 'bacgenome.report.homofilter'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        params = {
            'task_id': data.task_id,
            'homology_id': data.homology_id,
            'filter_type': data.filter_type,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'all_samples': data.all_samples
        }
        if hasattr(data, "filter"):
            params["filter"] = data.filter
        if hasattr(data,"filter_f"):
            params["filter_f"] = data.filter_f
        main_table_name = 'homofilter_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params, sort_keys=True, separators=(',', ':')))
        ]
        stat_file = self.bacgenome.get_orthmcl_file(data.task_id, data.homology_id)
        options = {
                    'specimens': data.all_samples,
                    'stat_file': stat_file,
                    'task_id': data.task_id,
                    'filter_type': data.filter_type
                   }

        if data.filter_type == 'core':
            options["select"] = data.all_samples
        elif data.filter_type == "unique":
            options["select"] = data.filter
        else:
            options["filter"] = data.filter_f
            options["select"] = data.filter

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_id = self.bacgenome.insert_main_table('compare_homofilter', mongo_data)
        update_info = {str(main_table_id): "compare_homofilter"}
        options["update_info"] = json.dumps(update_info)
        options["main_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(HomofilterAction, self).POST()
        # return spe_list
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
