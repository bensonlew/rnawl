# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.models.mongo.group_stat import GroupStat as G
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from bson import SON


class TwoGroupAction(MetaasvController):
    """
    Metaasv 两组分析
    """
    def __init__(self):
        super(TwoGroupAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'metaasv.report.two_group'
        task_type = 'workflow'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        # task_info = self.metaasv.get_task_info(otu_info['task_id'])
        project_sn = otu_info['project_sn']
        task_id = otu_info['task_id']
        main_table_name = 'DiffStatTwoGroup_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        groupname = json.loads(data.group_detail).keys()
        groupname.sort()
        category_name = ','.join(groupname)
        my_param = dict()
        my_param['asv_id'] = data.asv_id
        my_param['level_id'] = int(data.level_id)
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        if hasattr(data,'ci_test'):
            my_param['ci_test'] = data.ci_test
        my_param['multi_test'] = data.multi_test
        my_param['tail_test'] = data.tail_test
        my_param['test_method'] = data.test_method
        # my_param['methor'] = data.methor
        if hasattr(data,'coverage'):
            my_param['coverage'] = float(data.coverage)
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = '2'
        if hasattr(data, 'pair_id'):
            my_param['pair_id'] = data.pair_id
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        mongo_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "asv_id": data.asv_id if isinstance(data.asv_id, ObjectId) else ObjectId(data.asv_id),
            "group_id": data.group_id,
            "name": main_table_name,
            "level_id": int(data.level_id),
            "params": params,
            "status": "start",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "category_name": category_name
        }
        main_table_id = self.metaasv.insert_none_table('two_group')
        update_info = {str(main_table_id): 'two_group'}
        options = {
            "otu_file": data.asv_id,
            "params": params,
            "level": int(data.level_id),
            "test": data.test_method,
            "group_file": data.group_id,
            "correction": data.multi_test,
            "type": data.tail_test,
            "group_name": self.metaasv.get_group_name(data.group_id),
            "coverage": float(data.ci_test),
            "group_detail": data.group_detail,
            "category_name": category_name,
            "update_info": json.dumps(update_info),
            "main_id": str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        to_file = ["metaasv.export_otu_table_by_level(otu_file)", "metaasv.export_group_table_by_detail(group_file)"]

        if hasattr(data, 'pair_id'):  # 添加秩和检验配对样本文件
            options['signal_pair_file'] = data.pair_id
            options['signal_pair_id'] = str(data.pair_id)
            to_file.append("metaasv.export_signal_pair_sample(signal_pair_file)")

        self.set_sheet_data(name=task_name, options=options, main_table_name="DiffStatTwoGroup/" + main_table_name, module_type=task_type, to_file=to_file)
        task_info = super(TwoGroupAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传进来的参数是否正确
        """
        params_name = ['asv_id', 'level_id', 'group_detail', 'group_id',  'multi_test', 'tail_test', 'test_method', 'ci_test', 'submit_location']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数!")
        if int(data.level_id) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            success.append("level_id不在范围内")
        if data.multi_test not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
            success.append("多重检验方法不在范围内")
        if data.tail_test not in ["two.side", "greater", "less"]:
            success.append("检验类型不在范围内")
        if data.test_method not in ["chi", "fisher", "kru_H", "mann", "anova", "student", "welch", "signal"]:
            success.append("所选的分析检验方法不在范围内")
        if float(data.ci_test) not in [0.90, 0.95, 0.98, 0.99, 0.999]:
            success.append('置信区间的置信度coverage不在范围值内')
        table_dict = json.loads(data.group_detail)
        if len(table_dict) != 2 or data.group_id == 'all':
            success.append("分析只适用于分组方案的分组类别数量为2的情况！")
        if not isinstance(table_dict, dict):
            success.append("传入的table_dict不是一个字典")
        return success
