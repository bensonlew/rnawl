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


class PcaAction(BacComparativeController):
    def __init__(self):
        super(PcaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['task_id', 'function', 'submit_location', 'task_type', 'group_detail', 'group_id', 'file_path']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.beta_diversity'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)
        if data.function not in ['cog', 'function', 'ko', 'pathway']:
            info = {'success': False, 'info': '%s参数缺少!' % data.function}
            return json.dumps(info)

        params = {
            'function': data.function,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_detail': group_detail_sort(data.group_detail),
            'group_id': data.group_id,
            'task_id': data.task_id
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Pca_' + data.function + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
        main_id = self.bac_comparative.insert_main_table('pca', mongo_data)
        update_info[str(main_id)] = 'pca'
        options = {'function': data.function,
                   'update_info': json.dumps(update_info),
                   # 'params': params,
                   'main_id': str(main_id),
                   # 'main_table_data': SON(mongo_data),
                   'group': data.group_id,
                   'group_detail': data.group_detail,
                   "analysis": 'pca',
                   'file_path': data.file_path
                   }
        to_file = "bac_comparative.export_group_table_by_detail(group)"
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="PCA/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            to_file=to_file)
        task_info = super(PcaAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)