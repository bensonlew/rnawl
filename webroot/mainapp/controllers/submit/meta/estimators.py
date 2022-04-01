# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
import datetime
from bson import SON
from collections import OrderedDict


class EstimatorsAction(MetaController):
    """

    """
    ESTIMATORS = ['ace', 'bergerparker', 'boneh', 'bootstrap', 'bstick', 'chao', 'coverage', 'default', 'efron','geometric', 'goodscoverage', 'heip', 'invsimpson', 'jack', 'logseries', 'npshannon', 'nseqs','qstat', 'shannon', 'shannoneven', 'shen', 'simpson', 'simpsoneven', 'smithwilson', 'sobs', 'solow','pd']

    def __init__(self):
        super(EstimatorsAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['otu_id', 'level_id', 'index_type', 'submit_location', "group_id"]
        index_types = []
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu], "code" : "C2201303"}
                return json.dumps(info)
        for index in data.index_type.split(','):
            index_types.append(index)
            if index not in self.ESTIMATORS:
                variables = []
                variables.append(index)
                info = {"success": False, "info": "指数类型不正确{}".format(index), 'code':'C2201301', 'variables':variables}
                return json.dumps(info)
        sort_index = data.index_type.split(',')
        # sort_index.sort()
        sort_index = ','.join(sort_index)

        task_name = 'meta.report.estimators'
        task_type = 'workflow'

        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'index_type': sort_index,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "group_id": data.group_id,
            "group_detail": group_detail_sort(data.group_detail)
        }

        # by houshuang 20191015 用于Beta多样性分析“按Alpha多样性指数着色”判断是否做过分析
        pca_check = {
            'group_detail': params_json["group_detail"],
            'level_id': params_json["level_id"],
            'otu_id': params_json["otu_id"]
        }
        # <<<

        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2201302'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])

        level_name = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"]

        main_table_name = 'Estimators' + level_name[int(data.level_id) - 1] + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            # by houshuang 20191015 用于Beta多样性分析“按Alpha多样性指数着色”判断是否做过分析
            ('index_type', sort_index),
            ('pca_check', json.dumps(pca_check, sort_keys=True, separators=(',', ':')))
            # <<<
        ]
        main_table_id = self.meta.insert_none_table('sg_alpha_diversity')
        update_info = {str(main_table_id): 'sg_alpha_diversity'}

        options = {
            "otu_file": data.otu_id,
            "otu_id": data.otu_id,
            "indices": data.index_type,
            "level": data.level_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "group_detail": data.group_detail,
            "group_id": data.group_id,
            'update_info': json.dumps(update_info),
            "est_id": str(main_table_id),
            'main_table_data': SON(mongo_data)
        }
        # self.to_file = 'meta.export_otu_table_by_level(otu_file)'
        to_file = 'meta.export_otu_table_by_detail(otu_file)'
        self.set_sheet_data(name=task_name, options=options, main_table_name="Estimators/" + main_table_name,
                            module_type=task_type, to_file=to_file)  # modified by hongdongxuan 20170322 在main_table_name前面加上文件输出的文件夹名
        task_info = super(EstimatorsAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
        # print self.returnInfo
        # return_info = json.loads(self.returnInfo)
        # return_info['content']["ids"]["index_types"] = index_types
        # print(return_info['content']["ids"]["index_types"])
        # print(return_info)
        # return json.dumps(return_info)
