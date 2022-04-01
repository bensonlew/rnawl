# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from mainapp.models.mongo.meta import Meta
from mainapp.models.mongo.group_stat import GroupStat as G


class RocNewAction(MetaController):
    def __init__(self):
        super(RocNewAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        otu_info = Meta().get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2203401'}
            return json.dumps(info)
        task_name = 'meta.report.roc_new'
        task_type = 'workflow'
        task_info = self.meta.get_task_info(otu_info['task_id'])
        main_table_name = 'RocNew_' + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        params_json = {
            'otu_id': data.otu_id,
            'level_id': data.level_id,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'second_group_id': data.second_group_id,
            'second_group_detail': group_detail_sort(data.second_group_detail),
            'submit_location': data.submit_location,
            'task_type': 'reportTask',
            'lefse_analysis': data.lefse_analysis,
            'two_group_analysis': data.two_group_analysis,
            'ran_for_analysis': data.ran_for_analysis,
            'roc_calc_method': data.roc_calc_method
        }
        if hasattr(data, 'env_id') and hasattr(data, 'env_labs'):
            params_json['env_id'] = data.env_id
            params_json['env_labs'] = data.env_labs
        if hasattr(data, 'lda_filter'):
                # and data.strict and data.start_level and data.end_level and data.lefse_cho and data.lefse_num:
            params_json['lda_filter'] = data.lda_filter
            params_json['strict'] = data.strict
            params_json['start_level'] = data.start_level
            params_json['end_level'] = data.end_level
            params_json['lefse_cho'] = data.lefse_cho
            params_json['lefse_num'] = data.lefse_num
        if hasattr(data, 'method_cal'):
            # and data.ci and data.q_test and data.q_test and data.two_group_cho and data.two_group_num:
            params_json['method_cal'] = data.method_cal
            params_json['ci'] = data.ci
            params_json['q_test'] = data.q_test
            params_json['two_group_cho'] = data.two_group_cho
            params_json['two_group_num'] = data.two_group_num
        if hasattr(data, 'tree_number'):
                # and data.Ran_for_num:
            params_json['tree_number'] = data.tree_number
            params_json['Ran_for_num'] = data.Ran_for_num
        if hasattr(data, 'roc_method_1'):
            # data.roc_method_1 and data.roc_method_2:
            params_json['roc_method_1'] = data.roc_method_1
            params_json['roc_method_2'] = data.roc_method_2
        if hasattr(data,'intersection'):
            params_json['intersection'] = data.intersection
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),  # maybe data.otu_id
            ('name', main_table_name),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('status', 'start'),
            ('desc', '个性化roc分析'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_table_id = self.meta.insert_main_table('sg_roc_new', mongo_data)
        update_info = {str(main_table_id): 'sg_roc_new'}
        options = {
            "lefse_analysis": data.lefse_analysis,
            "two_group_analysis": data.two_group_analysis,
            "ran_for_analysis": data.ran_for_analysis,
            "otu_id": data.otu_id,
            "lefse_otu": data.otu_id,
            "lefse_group": data.group_id,
            "group_detail": data.group_detail,
            "second_group_detail": data.second_group_detail,
            "group_name": G().get_group_name(data.group_id, lefse=True, second_group=data.second_group_detail),
            "two_ran_otu": data.otu_id,
            "two_ran_group": data.group_id,
            "level": data.level_id,
            "roc_calc_method": data.roc_calc_method,
            "main_id": str(main_table_id),
            "update_info": json.dumps(update_info),
        }
        to_file = ["meta.export_otu_table(lefse_otu)", "meta.export_cascading_table_by_detail(lefse_group)",
                   "meta.export_otu_table_by_detail(two_ran_otu)", 'meta.export_group_table_by_detail(two_ran_group)']
        # if data.env_id and data.env_labs:
        if hasattr(data, 'env_id') and hasattr(data, 'env_labs'):
            options['env_id'] = data.env_id
            options['env_labs'] = data.env_labs
            options['env_table'] = data.env_id
            to_file.append('env.export_env_table(env_file)')
        # if data.lda_filter and data.strict and data.start_level and data.end_level and data.lefse_cho and data.lefse_num:
        if hasattr(data, 'lda_filter'):
            options['lda_filter'] = data.lda_filter
            options['strict'] = data.strict
            options['start_level'] = data.start_level
            options['end_level'] = data.end_level
            options['lefse_cho'] = data.lefse_cho
            options['lefse_num'] = data.lefse_num
        # if data.method_cal and data.ci and data.q_test and data.q_test and data.two_group_cho and data.two_group_num:
        if hasattr(data, 'method_cal'):
            options['method_cal'] = data.method_cal
            options['ci'] = data.ci
            options['q_test'] = data.q_test
            options['two_group_cho'] = data.two_group_cho
            options['two_group_num'] = data.two_group_num
        if hasattr(data, 'tree_number'):
        # if data.tree_number and data.Ran_for_num:
            options['tree_number'] = data.tree_number
            options['Ran_for_num'] = data.Ran_for_num
        # if data.roc_method_1 and data.roc_method_2:
        if hasattr(data, 'roc_method_1'):
            options['roc_method_1'] = data.roc_method_1
            options['roc_method_2'] = data.roc_method_2
        if hasattr(data, 'intersection'):
            options['intersection'] = data.intersection
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(RocNewAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端的参数
        :param data:
        :return: success
        """
        params_name = ['lefse_analysis', 'two_group_analysis', 'ran_for_analysis', 'otu_id', 'group_id', 'group_detail',
                       'second_group_id', 'second_group_detail', 'level_id', 'submit_location', 'task_type',
                       'roc_calc_method']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数!")
        if hasattr(data, 'strict'):
            if int(data.strict) not in [1, 0]:
                info = {"success": False, "info": "严格性比较策略不在范围内！"}
                return json.dumps(info)
        if hasattr(data, 'lda_filter'):
            if float(data.lda_filter) > 4.0 or float(data.lda_filter) < -4.0:
                success.append("LDA阈值不在范围内")
        if hasattr(data, 'start_level'):
            if int(data.start_level) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
                success.append('起始分类水平不在范围内')
        if hasattr(data, 'start_level'):
            if int(data.end_level) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
                success.append('结束分类水平不在范围内')
        group_detail = json.loads(data.group_detail)
        if not isinstance(group_detail, dict):
            success.append("传入的group_detail不是一个字典")
        elif len(group_detail) != 2:
            success.append("一级分组方案必须选择两个分组!")
        if data.second_group_detail != '':
            second_group_detail = json.loads(data.second_group_detail)
            first = 0
            second = 0
            for i in group_detail.values():
                first += len(i)
            for n in second_group_detail.values():
                second += len(n)
            if not isinstance(second_group_detail, dict):
                success.append("传入的second_group_detail不是一个字典")
            if first != second:
                success.append("二级分组与一级分组的样本数不相同，请检查！")
        if hasattr(data, "ci"):
            if float(data.ci) < 0.0 or float(data.ci) > 1.0:
                success.append("显著水平不在范围内")
        if data.roc_calc_method == "MI":
            if not (hasattr(data, "roc_method_1")) and not (hasattr(data, "roc_method_2")):
                success.append("MI分析缺少参数roc_method!")
        return success