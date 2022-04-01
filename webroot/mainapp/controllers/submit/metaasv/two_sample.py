# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.models.mongo.group_stat import GroupStat as G
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from bson import SON


class TwoSampleAction(MetaasvController):
    """
    Metaasv 两样本比较
    """
    def __init__(self):
        super(TwoSampleAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'metaasv.report.two_sample'
        task_type = 'workflow'  # 可以不配置
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        # task_info = self.metaasv.get_task_info(otu_info['task_id'])
        project_sn = otu_info['project_sn']
        task_id = otu_info['task_id']
        main_table_name = 'DiffStatTwoSample_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        my_param = dict()
        my_param['asv_id'] = data.asv_id
        my_param['level_id'] = int(data.level_id)
        my_param['sample1'] = data.sample1
        my_param['sample2'] = data.sample2
        if hasattr(data,'ci_test'):
            my_param['ci_test'] = data.ci_test
        my_param['multi_test'] = data.multi_test
        my_param['tail_test'] = data.tail_test
        my_param['test_method'] = data.test_method
        if hasattr(data,'coverage'):
            my_param['methor'] = data.methor
        if hasattr(data,'coverage'):
            my_param['coverage'] = float(data.coverage)
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = '2'
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        mongo_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "asv_id": data.asv_id if isinstance(data.asv_id, ObjectId) else ObjectId(data.asv_id),
            "name": main_table_name,
            "level_id": int(data.level_id),
            "params": params,
            "status": "start",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        main_table_id = self.metaasv.insert_none_table('two_sample')
        update_info = {str(main_table_id): 'two_sample'}
        options = {
            "otu_file": data.asv_id,
            "level": data.level_id,
            "test": data.test_method,
            "correction": data.multi_test,
            "type": data.tail_test,
            "sample1": data.sample1,
            "sample2": data.sample2,
            "methor": data.methor,
            "coverage": float(data.coverage),
            "params": params,
            "update_info": json.dumps(update_info),
            "main_id": str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        to_file = "metaasv.export_otu_table_by_level(otu_file)"
        self.set_sheet_data(name=task_name, options=options, main_table_name="DiffStatTwoSample/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(TwoSampleAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传进来的参数是否正确
        """
        params_name = ['asv_id', 'level_id', 'sample1', 'sample2', 'multi_test', 'tail_test', 'test_method', 'methor',
                       'coverage', 'submit_location']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数!")
        if int(data.level_id) not in [3, 4, 5, 6, 7, 8, 9]:
            success.append("level_id不在范围内")
        if data.multi_test not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
            success.append("多重检验方法不在范围内")
        if data.tail_test not in ["two.side", "greater", "less"]:
            success.append("检验类型不在范围内")
        if data.test_method not in ["chi", "fisher"]:
            success.append("所选的分析检验方法不在范围内")
        if float(data.coverage) not in [0.90, 0.95, 0.98, 0.99, 0.999]:
            success.append('置信区间的置信度不在范围值内')
        if data.methor not in ["DiffBetweenPropAsymptoticCC", "DiffBetweenPropAsymptotic", "NewcombeWilson"]:
            success.append('计算置信区间的方法不在范围值内')
        sample_name = self.metaasv.get_otu_sample_name(data.asv_id)
        if data.sample1 not in sample_name or data.sample2 not in sample_name:
            success.append('所输入的样本名不在otu表里，请检查样本名')
        return success
