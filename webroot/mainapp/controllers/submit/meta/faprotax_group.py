# -*- coding: utf-8 -*-
# __author__ = 'zzg'

import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.models.mongo.group_stat import GroupStat as G
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from bson import SON


class FaprotaxGroupAction(MetaController):

    def __init__(self):
        super(FaprotaxGroupAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result), "code" : "C2203702"}
            return json.dumps(info)
        task_name = 'meta.report.faprotax_group'
        task_type = 'workflow'  # 可以不配置
        print data.faprotax_id
        faprotax_info = self.meta.get_faprotax_info(data.faprotax_id)
        print faprotax_info
        if not faprotax_info:
           info = {"success": False, "info": "faprotax分析不存在，请确认参数是否正确！!"}
           return json.dumps(info)
        table_dict = json.loads(data.group_detail)
        if data.group_method != "two":
            if len(table_dict) <= 2 or data.group_id == 'all':
                info = {"success": False, "info": "分析只适用于分组方案的分组类别数量大于等于3的情况！", 'code': 'C2202301'}
                return json.dumps(info)
        task_info = self.meta.get_task_info(faprotax_info['task_id'])
        main_table_name = 'faprotaxGroup_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        groupname = json.loads(data.group_detail).keys()
        groupname.sort()
        category_name = ','.join(groupname)
        my_param = dict()
        my_param['faprotax_id'] = data.faprotax_id
        my_param['otu_id'] = data.otu_id
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        if hasattr(data,'ci'):
            my_param['ci'] = float(data.ci)
        my_param['correction'] = data.correction
        if hasattr(data,'type'):
            my_param['type'] = data.type
        my_param['test'] = data.test
        my_param['methor'] = data.methor
        my_param['coverage'] = float(data.coverage)
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = int(data.task_type)
        my_param['group_method'] = data.group_method
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        mongo_data = {
            "type": data.group_method,
            "project_sn": task_info['project_sn'],
            "task_id": task_info['task_id'],
            "faprotax_id": ObjectId(data.faprotax_id),
            "otu_id": data.otu_id,
            "group_id": data.group_id,
            "name": main_table_name,
            "params": params,
            "status": "start",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "category_name": category_name
        }
        main_table_id = self.meta.insert_none_table('sg_faprotax_group')
        update_info = {str(main_table_id): 'sg_faprotax_group'}
        if data.group_method == "two":
            options = {
                "faprotax_table": data.faprotax_id,
                "faprotax_id": data.faprotax_id,
                "params": params,
                "group_method": data.group_method,
                "test": data.test,
                "group_file": data.group_id,
                "group_id": data.group_id,
                "correction": data.correction,
                # "ci": float(data.ci),
                "type": data.type,
                "group_name": G().get_group_name(data.group_id),
                "coverage": data.coverage,
                "group_detail": data.group_detail,
                "category_name": category_name,
                "update_info": json.dumps(update_info),
                "main_id": str(main_table_id),
                "main_table_data": SON(mongo_data)
            }
        else:
            options = {
                "faprotax_table": data.faprotax_id,
                "faprotax_id": data.faprotax_id,
                "params": params,
                "group_method": data.group_method,
                "group_file": data.group_id,
                "group_id": data.group_id,
                "correction": data.correction,
                "test": data.test,
                "group_name": G().get_group_name(data.group_id),
                "methor": data.methor,
                "coverage": data.coverage,
                "group_detail": data.group_detail,
                "category_name": category_name,
                "update_info": json.dumps(update_info),
                "main_id": str(main_table_id),
                'main_table_data': SON(mongo_data)
            }
        to_file = ["meta.export_faprotax_table_by_faprotax_id(faprotax_table)", "meta.export_group_table_by_detail(group_file)"]

        self.set_sheet_data(name=task_name, options=options, main_table_name="faprotaxGroup/" + main_table_name, module_type=task_type, to_file=to_file)
        task_info = super(FaprotaxGroupAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传进来的参数是否正确
        """
        params_name = ['faprotax_id', 'group_detail', 'group_id',  'correction', 'test', 'coverage', 'submit_location']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数!")
        if data.correction not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
            success.append("多重检验方法不在范围内")
        if hasattr(data, 'type'):
            if data.type not in ["two.side", "greater", "less"]:
                success.append("检验类型不在范围内")
        if data.test not in ["chi", "fisher", "kru_H", "mann", "anova", "student", "welch", "signal","kru_H", "anova"]:
            success.append("所选的分析检验方法不在范围内")
        if float(data.coverage) not in [0.90, 0.95, 0.98, 0.99, 0.999]:
            success.append('置信区间的置信度coverage不在范围值内')
        table_dict = json.loads(data.group_detail)
        if data.group_method == "two":
            if len(table_dict) != 2 or data.group_id == 'all':     #modified by hongdongxuan 20170313
                success.append("选择的分析只适用于分组方案的分组类别数量为2的情况！")
        if not isinstance(table_dict, dict):
            success.append("传入的table_dict不是一个字典")
        return success
