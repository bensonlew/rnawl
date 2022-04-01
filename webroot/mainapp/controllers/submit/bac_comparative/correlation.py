# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20191122
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from bson import SON


class CorrelationAction(BacComparativeController):
    def __init__(self):
        super(CorrelationAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id', 'function', 'submit_location', 'task_type', 'group_detail', 'group_id', 'file_path', "dis_method", "hcluster_method", "correlation", "is_group"]
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

        if not data.dis_method:
            info = {'success': False, 'info': '%s参数缺少!' % data.dis_method}
            return json.dumps(info)

        if not data.hcluster_method in ['complete', 'average', 'single']:
            info = {'success': False, 'info': '%s参数不在["complete", "average", "single"]!' % data.hcluster_method}
            return json.dumps(info)

        if not data.is_group in ['yes', 'no']:
            info = {'success': False, 'info': '%s参数不在["yes","no"]中!' % data.is_group}
            return json.dumps(info)

        if not data.correlation in ['pearson', 'spearman', 'kendalltau']:
            info = {'success': False, 'info': '%s参数不在["pearson", "spearman", "kendalltau"]中!' % data.correlation}
            return json.dumps(info)
        correlation = ''
        if data.correlation in ['pearson']:
            correlation = 'pearsonr'
        elif data.correlation in ['spearman']:
            correlation = 'spearmanr'

        params = {
            'function': data.function,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_detail': group_detail_sort(data.group_detail),
            'group_id': data.group_id,
            "dis_method": data.dis_method,
            "hcluster_method": data.hcluster_method,
            "correlation": data.correlation,
            'task_id': data.task_id,
            'is_group': data.is_group,
        }
        if hasattr(data, "group_method"):
            params["group_method"] = data.group_method
            if not data.group_method in ["",'sum', 'average', 'middle']:
                info = {'success': False, 'info': '%s参数不在["sum", "average", "middle", ""]!' % data.group_method}
                return json.dumps(info)

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Correlation_' + data.function + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
        main_id = self.bac_comparative.insert_main_table('correlation', mongo_data)
        update_info[str(main_id)] = 'correlation'
        options = {'function': data.function,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   'group': data.group_id,
                   'group_detail': data.group_detail,
                   "analysis": 'correlation',
                   'file_path': data.file_path,
                   "dis_method": data.dis_method,
                   "hcluster_method": data.hcluster_method,
                   "corr_method": correlation
                   }
        if hasattr(data, "group_method"):
            options["group_method"] = data.group_method

        if hasattr(data, "is_group"):
            options["is_group"] = data.is_group

        to_file = "bac_comparative.export_group_table_by_detail(group)"
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Correlation/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            to_file=to_file)
        task_info = super(CorrelationAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)