# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
import datetime
from bson import SON


class AlphaDiversityDiffAction(MetaasvController):
    """
    Alpha 多样性指数组间差异检验
    """
    def __init__(self):
        super(AlphaDiversityDiffAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['alpha_diversity_id', 'group_detail', 'group_id', 'submit_location', 'test_method']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu]}
                return json.dumps(info)
        try:
            group = json.loads(data.group_detail)
        except ValueError:
            variables = [data.group_detail]
            info = {'success': False, 'info': 'group_detail格式不正确!:%s' % data.group_detail, 'variables': variables}
            return json.dumps(info)
        if data.test_method == "signal":
            temp = [len(i) for i in group.values()]
            for i in range(len(temp)-1):
                for j in range(i+1, len(temp)):
                    if temp[i] != temp[j]:
                        info = {'success': False, 'info': '存在两分组样本数不相等'}
                        return json.dumps(info)
        table_dict = json.loads(data.group_detail)
        if not isinstance(table_dict, dict):
            info = {"success": False, "info": "传入的group_detail不是字典"}
            return json.dumps(info)
        if len(table_dict) < 2:
            info = {"success": False, "info": "请选择至少两组以上的分组"}
            return json.dumps(info)
        for key in table_dict:
            if len(table_dict[key]) < 2:
                info = {"success": False, "info": "每组至少有两个样本"}
                return json.dumps(info)

        task_name = 'metaasv.report.alpha_compare'
        task_type = 'workflow'

        params_json = {
            'alpha_diversity_id': data.alpha_diversity_id,
            # "otu_id": data.otu_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "group_id": data.group_id,
            "group_detail": group_detail_sort(data.group_detail),
            "analysis_type": data.analysis_type,
        }
        if hasattr(data, "test_method"):
            params_json["test_method"] = data.test_method

        if hasattr(data, "multi_test"):
            params_json["multi_test"] = data.multi_test

        if hasattr(data, "post_hoc"):
            params_json["post_hoc"] = data.post_hoc

        if hasattr(data, "coverage"):
            params_json["coverage"] = data.coverage

        if hasattr(data, "ci_test"):
            params_json["ci_test"] = data.ci_test
        if hasattr(data, "tail_test"):
            params_json["tail_test"] = data.tail_test

        otu_info = self.metaasv.get_diversity_table_info(data.alpha_diversity_id)
        asv_id = otu_info["asv_id"]
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        main_table_name = 'EstimatorStat_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', asv_id),
            ('alpha_diversity_id', ObjectId(data.alpha_diversity_id)),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_none_table('alpha_diversity_diff')
        update_info = {str(main_table_id): 'alpha_diversity_diff'}
        options = {
            "est_table": data.alpha_diversity_id,
            "asv_id": str(asv_id),
            "group_table": data.group_id,
            'update_info': json.dumps(update_info),
            "est_id": data.alpha_diversity_id,
            "group_detail": data.group_detail,
            "main_id": str(main_table_id),
            "group_name": self.metaasv.get_group_name(data.group_id),
            'main_table_data': SON(mongo_data)
        }
        if hasattr(data, "test_method"):
            options["test"] = data.test_method
        if hasattr(data, "multi_test"):
            options["correction"] = data.multi_test
        if hasattr(data, "post_hoc"):
            options["methor"] = data.post_hoc
        if hasattr(data, "coverage"):
            options["coverage"] = float(data.coverage)
        if hasattr(data, "ci_test"):
            options["ci"] = float(data.ci_test)
        if hasattr(data, "tail_test"):
            options["type"] = data.tail_test
        if hasattr(data, "analysis_type"):
            options["analysis"] = data.analysis_type

        to_file = ["metaasv.export_est_table(est_table)", "metaasv.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="Estimators/" + main_table_name,
                            module_type=task_type, to_file=to_file, main_id=str(asv_id),collection_name="asv")
        task_info = super(AlphaDiversityDiffAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        if not task_info["success"]:
            task_info["info"] = "程序运行出错，请检查输入的多样性指数表是否存在异常（样本值完全相同或是存在NA值等情况）"
        return json.dumps(task_info)
