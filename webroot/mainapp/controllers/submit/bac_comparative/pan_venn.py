# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20191122
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from bson import SON
from mainapp.libs.param_pack import group_detail_sort


class PanVennAction(BacComparativeController):
    def __init__(self):
        super(PanVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'pan_id', 'submit_location', 'task_type','group_detail', 'group_id', 'cluster_path']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.pan_venn'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)

        params = {
            'pan_id': data.pan_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_detail' : group_detail_sort(data.group_detail),
            'group_id': data.group_id,
            'task_id': data.task_id
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Pan_Venn' + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.bac_comparative.insert_main_table('pan_venn', mongo_data)
        update_info[str(main_id)] = 'pan_venn'
        options = {'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   'group_table': data.group_id,
                   'group_detail': data.group_detail,
                   'cluster_file': data.cluster_path
                   }
        to_file = "bac_comparative.export_group_table_by_detail(group_table)"
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Pan_Venn/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            to_file=to_file)
        task_info = super(PanVennAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)