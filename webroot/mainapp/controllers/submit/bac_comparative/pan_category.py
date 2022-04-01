# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20191122
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.bac_comparative import BacComparative
from bson import SON


class PanCategoryAction(BacComparativeController):
    def __init__(self):
        super(PanCategoryAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # print(data)
        default_argu = ['task_id', 'pan_id','submit_location', 'task_type', 'cluster_path', "category_id"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.pan_category'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)

        category_id = data.category_id
        bac_genome = BacComparative()
        (category, percent, category_name) = bac_genome.get_pan_category(category_id)
        params = {
            'pan_id': data.pan_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'category_id': category_id,
            'task_id': data.task_id
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = "Pan_Genome_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            # ("category_scheme", category_scheme),
            ("pan_category", category_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.bac_comparative.insert_main_table('pan_genome', mongo_data)
        update_info[str(main_id)] = 'pan_genome'
        options = {'from_main_id': data.pan_id,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   'category': category,
                   "percent": percent,
                   "category_name": category_name,
                   'cluster_file': data.cluster_path,
                   # "class": category_scheme,
                   "pan_group_id": category_id,
                   'pan_category_names': self.bac_comparative.get_value_of_key('pan_group',
                                                                        category_id,
                                                                        'category_names')
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Pan_Genome/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(PanCategoryAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)