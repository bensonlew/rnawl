# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
import web
import json
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.param_pack import group_detail_sort
# import bson.errors.InvalidId
import types
from mainapp.models.mongo.meta import Meta
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig
from bson import SON


class MultiAnalysisAction(MetaController):
    def __init__(self):
        super(MultiAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['analysis_type', 'otu_id', 'level_id', 'submit_location', 'group_id', 'group_detail']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'PARAMETERS MISSING: %s' , "variables":[ argu], "code" : "C2200316"}
                return json.dumps(info)
        task_name = 'meta.report.beta_multi_analysis'
        task_type = 'workflow'
        meta = Meta()
        otu_info = meta.get_otu_table_info(data.otu_id)
        if hasattr(data, 'scale'):
            if data.scale not in ["T","F"]:
                info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of scale, T or F expected!', "code" : "C2200317"}
                return json.dumps(info)
        if not otu_info:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of otu_id, not found!", "code" : "C2200318"}
            return json.dumps(info)
        task_info = meta.get_task_info(otu_info['task_id'])
        group_detail = group_detail_sort(data.group_detail)
        main_table_name = MultiAnalysisAction.get_main_table_name(data.analysis_type)
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'analysis_type': data.analysis_type,
            'submit_location': data.submit_location,
            'task_type': data.task_type,
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
            ('otu_id', ObjectId(data.otu_id)),
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
            info = {'success': False, 'info': '样本数量少于3，不可进行此分析！', 'code':'C2200303'}
            return json.dumps(info)
        if data.analysis_type == 'pca':
            if hasattr(data, 'scale'):
                params_json['scale'] = data.scale
            if hasattr(data, 'env_id'):
                params_json['env_id'] = data.env_id
                env_id = self.check_objectid(data.env_id)
                if not env_id:
                    variables = []
                    variables.append(data.env_id)
                    info = {'success': False, 'info': 'env_id格式:%s不正确，无法转换为ObjectId格式！' % data.env_id, 'code':'C2200301', 'variables':variables}
                    return json.dumps(info)
                if hasattr(data, 'env_labs'):
                    params_json['env_labs'] = data.env_labs
                    env_labs = data.env_labs
                else:
                    info = {'success': False, 'info': '没有选择任何环境因子列', 'code':'C2200302'}
                    return json.dumps(info)
            # by houshuang 20190924 增加组间差异检验>>>
            if hasattr(data, 'diff_test_method') and hasattr(data, 'change_times'):
                if data.diff_test_method and data.diff_test_method != "none":
                    params_json['diff_test_method'] = data.diff_test_method
                params_json['change_times'] = str(data.change_times)
            # <<<
        elif data.analysis_type == 'pcoa' or data.analysis_type == 'nmds':
            if not hasattr(data, 'distance_algorithm'):
                info = {'success': False, 'info': 'distance_algorithm参数缺少!', 'code':'C2200304'}
                return json.dumps(info)
            params_json['distance_algorithm'] = data.distance_algorithm
            dist_method = data.distance_algorithm
            # by houshuang 20190924 增加组间差异检验>>>
            if hasattr(data, 'diff_test_method') and hasattr(data, 'change_times'):
                if data.diff_test_method and data.diff_test_method != "none":
                    params_json['diff_test_method'] = data.diff_test_method
                params_json['change_times'] = data.change_times
            # <<<
        elif data.analysis_type == 'dbrda':
            if not hasattr(data, 'distance_algorithm'):
                info = {'success': False, 'info': 'distance_algorithm参数缺少!', 'code':'C2200305'}
                return json.dumps(info)
            params_json['distance_algorithm'] = data.distance_algorithm
            dist_method = data.distance_algorithm
            if hasattr(data, 'env_id'):
                params_json['env_id'] = data.env_id
                env_id = self.check_objectid(data.env_id)
                if not env_id:
                    variables = []
                    variables.append(data.env_id)
                    info = {'success': False, 'info': 'env_id格式:%s不正确，无法转换为ObjectId格式！' % data.env_id, 'code':'C2200306', 'variables':variables}
                    return json.dumps(info)
                if hasattr(data, 'env_labs'):
                    params_json['env_labs'] = data.env_labs
                    env_labs = data.env_labs
                else:
                    info = {'success': False, 'info': '没有选择任何环境因子列', 'code':'C2200307'}
                    return json.dumps(info)
            else:
                # info = {'success': False, 'info': 'dbrda分析缺少参数:env_id!'}
                info = {'success': False, 'info': 'dbrda分析缺少环境因子参数!', 'code':'C2200308'}  #modified by hongdongxuan 20170310
                return json.dumps(info)
        elif data.analysis_type == 'rda_cca':
            if hasattr(data, 'env_id'):
                params_json['env_id'] = data.env_id
                env_id = self.check_objectid(data.env_id)
                if not env_id:
                    variables = []
                    variables.append(data.env_id)
                    info = {'success': False, 'info': 'env_id格式:%s不正确，无法转换为ObjectId格式！' % data.env_id, 'code':'C2200309', 'variables':variables}
                    return json.dumps(info)
                if hasattr(data, 'env_labs'):
                    params_json['env_labs'] = data.env_labs
                    env_labs = data.env_labs
                else:
                    info = {'success': False, 'info': '没有选择任何环境因子列', 'code':'C2200310'}
                    return json.dumps(info)
            else:
                # info = {'success': False, 'info': 'rda_cca分析缺少参数:env_id!'}
                info = {'success': False, 'info': 'rda_cca分析缺少环境因子参数!', 'code':'C2200311'}  #modified by hongdongxuan 20170310
                return json.dumps(info)
        elif data.analysis_type == 'plsda':
            try:
                group = json.loads(data.group_detail)
            except ValueError:
                variables = []
                variables.append(data.group_detail)
                info = {'success': False, 'info': 'group_detail格式不正确!:%s' % data.group_detail, 'code':'C2200312','variables':variables}
                return json.dumps(info)
            # params_json['group_detail'] = group_detail_sort(data.group_detail)
            if len(group) < 2:
                info = {'success': False, 'info': '分析只适用于分组方案的分组类别数量大于等于2的情况！', 'code':'C2200313'}
                return json.dumps(info)
            samples = reduce(lambda x, y: x + y, group.values())
            if len(samples) == len(set(samples)):
                pass
            else:
                info = {'success': False, 'info': '不同分组存在相同的样本id', 'code':'C2200314'}
                return json.dumps(info)
        else:
            variables = []
            variables.append(data.analysis_type)
            info = {'success': False, 'info': '不正确的分析方法:%s' % data.analysis_type, 'code':'C2200315', 'variables':variables}
            return json.dumps(info)
        options = {
            'analysis_type': data.analysis_type,
            'otu_file': data.otu_id,
            'otu_id': data.otu_id,
            'level': int(data.level_id),
            'dist_method': dist_method,
            'env_labs': env_labs,
            'group_id': data.group_id,
            'group_detail': data.group_detail,
            'good_group': good_group,
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
        }
        to_file = ['meta.export_otu_table_by_detail(otu_file)']
        mongo_data.append(('env_id', env_id))
        # by houshuang 20190924 增加组间差异检验>>>
        if hasattr(data, 'diff_test_method') and hasattr(data, 'change_times'):
            try:
                group = json.loads(data.group_detail)
            except ValueError:
                variables = []
                variables.append(data.group_detail)
                info = {'success': False, 'info': 'group_detail格式不正确!:%s' % data.group_detail, 'variables': variables, "code" : "C2200319"}
                return json.dumps(info)
            if len(group) < 2:
                info = {'success': False, 'info': '分析只适用于分组方案的分组类别数量大于等于2的情况！', "code" : "C2200320"}
                return json.dumps(info)
            for i in group.values():
                if len(i) < 3:
                    info = {'success': False, 'info': '分析只适用于每个分组内样本数均大于等于2的情况！', "code" : "C2200321"}
                    return json.dumps(info)
            try:
                int(data.change_times)
            except ValueError:
                variables = []
                variables.append(data.change_times)
                info = {'success': False, 'info': 'permutations格式应该为数字!:%s' % data.change_times, 'variables': variables, "code" : "C2200322"}
                return json.dumps(info)
            if not (9 < int(data.change_times) < 10001):
                variables = []
                variables.append(data.change_times)
                info = {'success': False, 'info': '置换次数应该在[10-10000]之间:%s' % data.change_times, 'variables': variables, "code" : "C2200323"}
                return json.dumps(info)
            if data.diff_test_method not in ['adonis', 'anosim']:
                variables = []
                variables.append(data.diff_test_method)
                info = {'success': False, 'info': '组间差异检验方法不存在:%s' % data.diff_test_method, 'variables': variables, "code" : "C2200324"}
                return json.dumps(info)
            options['diff_test_method'] = data.diff_test_method
            options['change_times'] = data.change_times
        # <<<
        if hasattr(data, 'scale'):
            options['scale'] = data.scale
        if env_id:
            mongo_data.append(('env_labs', data.env_labs))
            to_file.append('env.export_env_table(env_file)')
            options['env_file'] = data.env_id
            options['env_id'] = data.env_id
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        # main_table_id = meta.insert_main_table('sg_beta_multi_analysis', mongo_data)
        main_table_id = meta.insert_none_table('sg_beta_multi_analysis')
        update_info = {str(main_table_id): 'sg_beta_multi_analysis'}
        options['update_info'] = json.dumps(update_info)
        options['main_id'] = str(main_table_id)
        if data.analysis_type in ['plsda', 'pca', 'pcoa', 'nmds', 'rda_cca', 'dbrda']:  # by houshuang 20190924 增加dbrda和rda_cca
            to_file.append('meta.export_group_table_by_detail(group_file)')
            options['group_file'] = data.group_id
        options['main_table_data'] = SON(mongo_data)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(MultiAnalysisAction, self).POST()
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
