# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
import types
from bson import ObjectId
from bson.errors import InvalidId
from bson import SON
from mainapp.models.mongo.metaasv import Metaasv
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig


class PermanovaAction(MetaasvController):
    """
    Metaasv Permanova样本菌群分析
    """
    def __init__(self):
        super(PermanovaAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'level_id', 'distance_method','change_times', 'group_detail','task_type','submit_location']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        try:
            group = json.loads(data.group_detail)
        except ValueError:
            variables = []
            variables.append(data.group_detail)
            info = {'success': False, 'info': 'group_detail格式不正确!:%s' % data.group_detail}
            return json.dumps(info)
        try:
            int(data.change_times)
        except ValueError:
            variables = []
            variables.append(data.change_times)
            info = {'success': False, 'info': 'change_times格式应该为数字!:%s' % data.change_times}
            return json.dumps(info)
        if not (9 < int(data.change_times) < 10000):
            variables = []
            variables.append(data.change_times)
            info = {'success': False, 'info': '置换次数应该在[10-10000]之间:%s' % data.change_times}
            return json.dumps(info)
        if (not hasattr(data, 'env_id')) and (not hasattr(data, 'group_detail_list')):
            info = {'success': False, 'info': '必须提供分组文件或环境因子文件!'}
            return json.dumps(info)
        distance_name = ['euclidean', 'binary_euclidean', 'manhattan', 'binary_manhattan', 'gowerM', 'binary_gowerM',
                                      'altGower', 'binary_altGower', 'canberraNZ', 'binary_canberraNZ', 'bray_curtis', 'binary_bray_curtis',
                                      'kulczynski', 'binary_kulczynski', 'morisita_horn', 'binary_morisita_horn', 'morisita',
                                      'binomial', 'binary_binomial', 'cao', 'binary_cao', 'chao', 'jaccard', 'binary_jaccard',
                                      'raup_crick', 'mountford', 'mahalanobis']
        if data.distance_method not in distance_name:
            info = {'success': False, 'info': '不支持所选的距离算法!'}
            return json.dumps(info)

        task_name = 'metaasv.report.permanova'
        module_type = 'workflow'  # 可以不配置
        main_table_name = 'Permanova_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        meta = Metaasv()
        otu_info = meta.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = meta.get_task_info(otu_info['task_id'])
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'distance_method': data.distance_method,
            'change_times': data.change_times,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        if hasattr(data, 'group_detail_list'):
            params_json['group_detail_list'] = eval(data.group_detail_list)
        if hasattr(data, 'env_id'):
            params_json['env_id'] = data.env_id
            env_id = self.check_objectid(data.env_id)
            if not env_id:
                variables = []
                variables.append(data.env_id)
                info = {'success': False, 'info': 'env_id格式:%s不正确，无法转换为ObjectId格式！' % data.env_id}
                return json.dumps(info)
            if hasattr(data, 'env_labs'):
                params_json['env_labs'] = data.env_labs
            else:
                info = {'success': False, 'info': '没有选择任何环境因子列'}
                return json.dumps(info)
        to_file = []
        main_table_id = meta.insert_none_table('permanova')
        update_info = {str(main_table_id): 'permanova'}
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        options = {
            'otu_file': data.asv_id,
            'asv_id': data.asv_id,
            'level': int(data.level_id),
            'method': data.distance_method,
            # 'task_type': data.task_type,
            'group_detail': data.group_detail,
            'permutation': int(data.change_times),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        if hasattr(data, 'group_detail_list'):
            options['group_detail_list'] = data.group_detail_list
        if hasattr(data, 'env_id'):
            options['env_id'] = data.env_id
            if hasattr(data, 'env_labs'):
                options['env_labs'] = data.env_labs
        options['group_table'] = ''
        options['envtable'] = ''
        to_file.append('metaasv.export_otu_table_by_level(otu_file)')
        to_file.append('metaasv.export_sample_list_by_detail(group_table)')
        to_file.append('metaasv.export_env_group_table(envtable)')

        self.set_sheet_data(name=task_name, options=options, main_table_name="Permanova/" + main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(PermanovaAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def check_objectid(self, in_id):
        """检查一个id是否可以被ObjectId"""
        if isinstance(in_id, types.StringTypes):
            try:
                in_id = ObjectId(in_id)
            except InvalidId:
                return False
        elif isinstance(in_id, ObjectId):
            pass
        else:
            return False
        return in_id
