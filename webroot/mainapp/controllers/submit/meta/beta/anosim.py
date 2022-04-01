# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.models.mongo.meta import Meta
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig


class AnosimAction(MetaController):
    def __init__(self):
        super(AnosimAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['otu_id', 'level_id', 'distance_algorithm',
                        'permutations', 'group_id', 'group_detail',
                        'submit_location']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing: %s' % argu}
                return json.dumps(info)
        try:
            group = json.loads(data.group_detail)
        except ValueError:
            variables = []
            variables.append(data.group_detail)
            info = {'success': False, 'info': 'group_detail格式不正确!:%s' % data.group_detail, 'code':'C2200101', 'variables':variables}
            return json.dumps(info)
        try:
            int(data.permutations)
        except ValueError:
            variables = []
            variables.append(data.permutations)
            info = {'success': False, 'info': 'permutations格式应该为数字!:%s' % data.permutations, 'code':'C2200102', 'variables':variables}
            return json.dumps(info)
        if not (9 < int(data.permutations) < 10000):
            variables = []
            variables.append(data.permutations)
            info = {'success': False, 'info': '置换次数应该在[10-10000]之间:%s' % data.permutations, 'code':'C2200103', 'variables':variables}
            return json.dumps(info)
        if len(group) < 2:
            info = {'success': False, 'info': '分析只适用于分组方案的分组类别数量大于等于2的情况！', 'code':'C2200104'}
            return json.dumps(info)
        samples = reduce(lambda x, y: x + y, group.values())
        if len(samples) == len(set(samples)):
            pass
        else:
            info = {'success': False, 'info': '不同分组存在相同的样本id', 'code':'C2200105'}
            return json.dumps(info)
        if len(samples) <= len(group):
            info = {'success': False, 'info': '不可每个组都只含有一个样本', 'code':'C2200106'}
            return json.dumps(info)

        task_name = 'meta.report.anosim'
        task_type = 'workflow'  # 可以不配置
        main_table_name = 'Anosim_Adonis_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        meta = Meta()
        otu_info = meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2200107'}
            return json.dumps(info)
        task_info = meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'distance_algorithm': data.distance_algorithm,
            'permutations': data.permutations,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = meta.insert_none_table('sg_beta_multi_anosim')
        update_info = {str(main_table_id): 'sg_beta_multi_anosim'}
        options = {
            'otu_file': data.otu_id,
            'otu_id': data.otu_id,
            'level': int(data.level_id),
            'method': data.distance_algorithm,
            'group_file': data.group_id,
            'group_detail': data.group_detail,
            'permutations': int(data.permutations),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        to_file = ['meta.export_otu_table_by_level(otu_file)',
                   'meta.export_group_table_by_detail(group_file)']
        self.set_sheet_data(name=task_name, options=options, main_table_name="Anosim_Adonis/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(AnosimAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
