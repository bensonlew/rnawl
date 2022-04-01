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
from mainapp.libs.signature import check_sig
from bson import SON
from mainapp.controllers.project.metaasv_controller import MetaasvController



class BetaDiversityAction(MetaasvController):
    """
    Metaasv PCA、PCoA、NMDS分析
    """
    def __init__(self):
        super(BetaDiversityAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['analysis_type', 'asv_id', 'level_id', 'submit_location', 'group_id', 'group_detail']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'PARAMETERS MISSING: %s'}
                return json.dumps(info)
        task_name = 'metaasv.report.beta_diversity'
        task_type = 'workflow'
        meta = Metaasv()
        otu_info = meta.get_otu_table_info(data.asv_id)
        if hasattr(data, 'scale'):
            if data.scale not in ["T","F"]:
                info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of scale, T or F expected!'}
                return json.dumps(info)
        if not otu_info:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of otu_id, not found!"}
            return json.dumps(info)
        task_info = meta.get_task_info(otu_info['task_id'])
        group_detail = group_detail_sort(data.group_detail)
        main_table_name = BetaDiversityAction.get_main_table_name(data.analysis_type)
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'analysis_type': data.analysis_type,
            'submit_location': data.submit_location,
            'task_type': str(data.task_type),
            'group_id': data.group_id,
            'group_detail': group_detail
        }
        dist_method = ''
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
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
        if data.analysis_type == 'pca':
            if hasattr(data, 'scale'):
                params_json['scale'] = data.scale
            if hasattr(data, 'diff_test_method') and hasattr(data, 'change_times'):
                params_json['diff_test_method'] = data.diff_test_method
                params_json['change_times'] = data.change_times
        elif data.analysis_type == 'pcoa' or data.analysis_type == 'nmds':
            if not hasattr(data, 'distance_algorithm'):
                info = {'success': False, 'info': 'distance_algorithm参数缺少!'}
                return json.dumps(info)
            params_json['distance_algorithm'] = data.distance_algorithm
            dist_method = data.distance_algorithm
            if hasattr(data, 'diff_test_method') and hasattr(data, 'change_times'):
                params_json['diff_test_method'] = data.diff_test_method
                params_json['change_times'] = data.change_times

        else:
            variables = []
            variables.append(data.analysis_type)
            info = {'success': False, 'info': '不正确的分析方法:%s' % data.analysis_type}
            return json.dumps(info)
        options = {
            'analysis_type': data.analysis_type,
            'otu_file': data.asv_id,
            'asv_id': data.asv_id,
            'level': int(data.level_id),
            'dist_method': dist_method,
            'group_id': data.group_id,
            'group_detail': data.group_detail,
            'good_group': good_group,
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
        }
        to_file = ['metaasv.export_otu_table_by_detail(otu_file)']
        if hasattr(data, 'diff_test_method') and hasattr(data, 'change_times'):
            try:
                group = json.loads(data.group_detail)
            except ValueError:
                variables = []
                variables.append(data.group_detail)
                info = {'success': False, 'info': 'group_detail格式不正确!:%s' % data.group_detail}
                return json.dumps(info)
            if len(group) < 2:
                info = {'success': False, 'info': '分析只适用于分组方案的分组类别数量大于等于2的情况！'}
                return json.dumps(info)
            for i in group.values():
                if len(i) < 3:
                    info = {'success': False, 'info': '分析只适用于每个分组内样本数均大于等于2的情况！'}
                    return json.dumps(info)
            try:
                int(data.change_times)
            except ValueError:
                variables = []
                variables.append(data.change_times)
                info = {'success': False, 'info': 'permutations格式应该为数字!:%s' % data.change_times}
                return json.dumps(info)
            if not (9 < int(data.change_times) < 10001):
                variables = []
                variables.append(data.change_times)
                info = {'success': False, 'info': '置换次数应该在[10-10000]之间:%s' % data.change_times}
                return json.dumps(info)
            if data.diff_test_method not in ['adonis', 'anosim']:
                variables = []
                variables.append(data.diff_test_method)
                info = {'success': False, 'info': '组间差异检验方法不存在:%s' % data.diff_test_method}
                return json.dumps(info)
            options['diff_test_method'] = data.diff_test_method
            options['change_times'] = data.change_times
        if hasattr(data, 'scale'):
            options['scale'] = data.scale

        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = meta.insert_none_table(data.analysis_type)
        update_info = {str(main_table_id): data.analysis_type}
        options['update_info'] = json.dumps(update_info)
        options['main_id'] = str(main_table_id)
        to_file.append('metaasv.export_group_table_by_detail(group_file)')
        options['group_file'] = data.group_id
        options['main_table_data'] = SON(mongo_data)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(BetaDiversityAction, self).POST()
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

    @staticmethod
    def get_main_table_name(analysis_type):
        time_now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        if analysis_type == 'pca':
            return 'PCA_' + time_now
        elif analysis_type == 'pcoa':
            return 'PCoA_' + time_now
        elif analysis_type == 'nmds':
            return 'NMDS_' + time_now
        elif analysis_type == 'plsda':
            return 'PLS-DA_' + time_now
        elif analysis_type == 'dbrda':
            return 'db-RDA_' + time_now
        elif analysis_type == 'rda_cca':
            return 'RDACCA_' + time_now
        else:
            raise Exception('PARAMETERS ERROR: wrong value of analysis_type')
