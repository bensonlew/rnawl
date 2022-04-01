# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.param_pack import group_detail_sort
import types
from mainapp.models.mongo.metaasv import Metaasv
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig
from bson import SON


class BetaDiversityEnvAction(MetaasvController):
    """
    metaasv beta多样性分析（RDA_CCA分析和dbrda分析）
    """
    def __init__(self):
        super(BetaDiversityEnvAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['analysis_type', 'asv_id', 'level_id', 'submit_location', 'group_id', 'group_detail']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'PARAMETERS MISSING: %s' , "variables":[ argu]}
                return json.dumps(info)
        task_name = 'metaasv.report.beta_diversity_env'
        task_type = 'workflow'
        meta = Metaasv()
        otu_info = meta.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of otu_id, not found!"}
            return json.dumps(info)
        task_info = meta.get_task_info(otu_info['task_id'])
        group_detail = group_detail_sort(data.group_detail)
        main_table_name = BetaDiversityEnvAction.get_main_table_name(data.analysis_type)
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'analysis_type': data.analysis_type,
            'submit_location': data.submit_location,
            'task_type': str(data.task_type),
            'group_id': data.group_id,
            'group_detail': group_detail
        }
        env_id = None
        env_labs = ''
        dist_method = ''
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('table_type', data.analysis_type),
            ('status', 'start'),
            ('group_id', group_id),
            ('desc', 'processing'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id))
        ]
        sample_len = sum([len(i) for i in group_detail.values()])
        group_len = sum([len(j) for j in group_detail.keys()])
        good_group = "F"   ## if there exists at least one group that has more than 3 samples, add by liulinmeng
        if group_len > 1:
            for k in group_detail.keys():
                if len(group_detail[k]) > 2:
                    good_group = "T"
                    break
        if sample_len < 3:
            info = {'success': False, 'info': '样本数量少于3，不可进行此分析！'}
            return json.dumps(info)
        if data.analysis_type == 'dbrda':
            if not hasattr(data, 'distance_algorithm'):
                info = {'success': False, 'info': 'distance_algorithm参数缺少!'}
                return json.dumps(info)
            params_json['distance_algorithm'] = data.distance_algorithm
            dist_method = data.distance_algorithm
            if hasattr(data, 'env_id'):
                params_json['env_id'] = data.env_id
                env_id = self.check_objectid(data.env_id)
                if not env_id:
                    variables = []
                    variables.append(data.env_id)
                    info = {'success': False, 'info': 'env_id格式:%s不正确，无法转换为ObjectId格式！' % data.env_id, 'variables':variables}
                    return json.dumps(info)
                if hasattr(data, 'env_labs'):
                    params_json['env_labs'] = data.env_labs
                    env_labs = data.env_labs
                else:
                    info = {'success': False, 'info': '没有选择任何环境因子列'}
                    return json.dumps(info)
            else:
                info = {'success': False, 'info': 'dbrda分析缺少环境因子参数!'}
                return json.dumps(info)
        elif data.analysis_type == 'rda_cca':
            if hasattr(data, 'env_id'):
                params_json['env_id'] = data.env_id
                env_id = self.check_objectid(data.env_id)
                if not env_id:
                    variables = []
                    variables.append(data.env_id)
                    info = {'success': False, 'info': 'env_id格式:%s不正确，无法转换为ObjectId格式！' % data.env_id, 'variables':variables}
                    return json.dumps(info)
                if hasattr(data, 'env_labs'):
                    params_json['env_labs'] = data.env_labs
                    env_labs = data.env_labs
                else:
                    info = {'success': False, 'info': '没有选择任何环境因子列'}
                    return json.dumps(info)
            else:
                info = {'success': False, 'info': 'rda_cca分析缺少环境因子参数!'}
                return json.dumps(info)
        options = {
            'analysis_type': data.analysis_type,
            'otu_file': data.asv_id,
            'asv_id': data.asv_id,
            'level': int(data.level_id),
            'dist_method': dist_method,
            'env_labs': env_labs,
            'group_id': data.group_id,
            'group_detail': data.group_detail,
            'group_file': data.group_id,
            'env_file': data.group_id,
            'good_group': good_group,
            # 'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
        }
        to_file = ['metaasv.export_otu_table_by_detail(otu_file)']
        mongo_data.append(('env_id', env_id))

        if env_id:
            mongo_data.append(('env_labs', data.env_labs))
            to_file.append('metaasv_env.export_env_table(env_file)')
            options['env_file'] = data.env_id
            options['env_id'] = data.env_id
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        if data.analysis_type in ['rda_cca']:
            main_table_id = meta.insert_none_table('rda_cca')
            update_info = {str(main_table_id): 'rda_cca'}
        if data.analysis_type in ['dbrda']:
            main_table_id = meta.insert_none_table('dbrda')
            update_info = {str(main_table_id): 'dbrda'}
        options['update_info'] = json.dumps(update_info)
        options['main_id'] = str(main_table_id)
        if data.analysis_type in ['rda_cca', 'dbrda']:
            to_file.append('metaasv.export_group_table_by_detail(group_file)')
            options['group_file'] = data.group_id
        options['main_table_data'] = SON(mongo_data)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(BetaDiversityEnvAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def check_objectid(self, in_id):
        """
        检查一个id是否可以被ObjectId
        """
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

    @staticmethod
    def get_main_table_name(analysis_type):
        time_now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        if analysis_type == 'dbrda':
            return 'dbRDA_' + time_now
        elif analysis_type == 'rda_cca':
            return 'RDACCA_' + time_now
        else:
            raise Exception('PARAMETERS ERROR: wrong value of analysis_type')
