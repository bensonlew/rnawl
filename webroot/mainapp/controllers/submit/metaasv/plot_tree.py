# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import SON


class PlotTreeAction(MetaasvController):
    def __init__(self):
        """
        Metaasv 个性化系统发生树分析
        """
        super(PlotTreeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['asv_id', 'level_id', 'color_level_id', 'group_id', 'group_detail', 'submit_location', 'topn', 'method']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu]}
                return json.dumps(info)
        if int(data.level_id) < int(data.color_level_id):
            info = {'success': False, 'info': '颜色设置水平必须高于进化树绘制水平!'}
            return json.dumps(info)
        if int(data.topn) < 2:
            if int(data.topn) != 0:
                info = {'success': False, 'info': '至少选择丰度高的物种2个及以上'}
                return json.dumps(info)
        task_name = 'metaasv.report.plot_tree'
        task_type = 'workflow'  # 可以不配置
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        task_id = otu_info['task_id']
        method = data.method   #'NJ'
        if method not in ['NJ', 'ML', 'MP']:
            info = {'success': False, 'info': "Value of  method parameters not in ['NJ', 'ML', 'MP']"}
            return json.dumps(info)

        if method in ['MP']:
            if int(data.topn) > 500:
                info = {'success': False, 'info': "top num can not more than 500"}
                return  json.dumps(info)
        elif method in ['NJ', 'ML']:
            if int(data.topn) > 500:
                info = {'success': False, 'info': "top num can not more than 500"}
                return  json.dumps(info)

        run_all_tree = 'start'

        params = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'color_level_id': int(data.color_level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            'topn': data.topn,
            'method' : method
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'PlotTree_' + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_id),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]

        if method == 'ML':
            settled_params = {"software": "IQ-TREE"}
            mongo_data.append(("settled_params", json.dumps(settled_params)))
        else:
            settled_params = {"software": "MEGA"}
            mongo_data.append(("settled_params", json.dumps(settled_params)))


        main_table_id = self.metaasv.insert_none_table('phylo_tree')
        update_info = {str(main_table_id): 'phylo_tree'}
        options = {
                   'asv_id': data.asv_id,
                   'level': int(data.level_id),
                   'color_level_id': int(data.color_level_id),
                   'group_id': data.group_id,
                   'update_info': json.dumps(update_info),
                   'group_detail': data.group_detail,
                   'params': params,
                   'main_id': str(main_table_id),
                   'topN': int(data.topn),
                   'main_table_data': SON(mongo_data),
                   'task_id' : task_info['task_id'],
                   'method' : method
                   }

        if run_all_tree == 'start':
            options["run_tree"] = 'part'
            options['otu_seq'] = data.asv_id
            to_file = ['metaasv.get_seq_by_select_samples_pick_pre_nums(otu_seq)']
        else:
            options['otu_table'] = data.asv_id
            options['sample_group'] = data.group_id
            options["run_tree"] = 'pick'
            to_file = ['metaasv.export_otu_table_by_detail(otu_table)', 'metaasv.export_group_table_by_detail(sample_group)']

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="PlotTree/" + main_table_name,
                            module_type=task_type,
                            to_file=to_file)
        task_info = super(PlotTreeAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
