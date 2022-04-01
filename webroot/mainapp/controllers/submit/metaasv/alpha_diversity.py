# -*- coding: utf-8 -*-
# __author__ = 'qingchen'
import web
import json
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
import datetime
from bson import SON


class AlphaDiversityAction(MetaasvController):
    """
    Alpha 多样性指数
    """
    ESTIMATORS = ['ace', 'bergerparker', 'boneh', 'bootstrap', 'bstick', 'chao', 'coverage', 'default', 'efron','geometric', 'goodscoverage', 'heip', 'invsimpson', 'jack', 'logseries', 'npshannon', 'nseqs','qstat', 'shannon', 'shannoneven', 'shen', 'simpson', 'simpsoneven', 'smithwilson', 'sobs', 'solow','pd']

    def __init__(self):
        super(AlphaDiversityAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'level_id', 'index_type', 'submit_location', "group_id"]
        index_types = []
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s'}
                return json.dumps(info)
        for index in data.index_type.split(','):
            index_types.append(index)
            if index not in self.ESTIMATORS:
                variables = []
                variables.append(index)
                info = {"success": False, "info": "指数类型不正确{}".format(index)}
                return json.dumps(info)
        sort_index = data.index_type.split(',')
        sort_index = ','.join(sort_index)

        task_name = 'metaasv.report.alpha_diversity'
        task_type = 'workflow'

        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'index_type': sort_index,
            "submit_location": data.submit_location,
            "task_type": str(data.task_type),
            "group_id": data.group_id,
            "group_detail": group_detail_sort(data.group_detail)
        }

        pca_check = {
            'group_detail': params_json["group_detail"],
            'level_id': params_json["level_id"],
            'asv_id': params_json["asv_id"]
        }

        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])

        level_name = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"]

        main_table_name = 'Estimators' + level_name[int(data.level_id) - 1] + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('index_type', sort_index),
            ('desc', '正在计算'),
            ('pca_check', json.dumps(pca_check, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_none_table('alpha_diversity')
        update_info = {str(main_table_id): 'alpha_diversity'}

        options = {
            "otu_file": data.asv_id,
            "asv_id": data.asv_id,
            "indices": data.index_type,
            "level": data.level_id,
            "task_type": data.task_type,
            "group_detail": data.group_detail,
            "group_table": data.group_id,
            'update_info': json.dumps(update_info),
            "main_id": str(main_table_id),
            'main_table_data': SON(mongo_data)
        }
        to_file = ['metaasv.export_otu_table_by_detail(otu_file)', "metaasv.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="Estimators/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(AlphaDiversityAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
