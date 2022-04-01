# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20191230

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig


class Picrust2Action(MetaController):
    def __init__(self):
        super(Picrust2Action, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['otu_id', 'submit_location', "group_id", "group_detail", "group_method", "task_type", "analysis_type", "hsp"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "parameters missing:%s" % argu}
                return json.dumps(info)
        if data.group_method not in ["", "sum", "average", "middle"]:
            info = {"success": False, "info": "对分组样本计算方式:%s错误!" % data.group_method}
            return json.dumps(info)
        task_name = 'meta.report.picrust2_predict'
        task_type = 'workflow'
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info["task_id"])
        params_json = {
            'otu_id': data.otu_id,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'group_method': data.group_method,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            "hsp": data.hsp,
            "analysis_type": data.analysis_type
        }
        main_table_name = "PICRUSt2_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', 'PICRUSt2功能预测'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('params', json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.meta.insert_none_table('sg_picrust2')
        update_info = {str(main_table_id): 'sg_picrust2'}
        options = {
            "update_info": json.dumps(update_info),
            "otu_table": data.otu_id,
            "otu_fasta": data.otu_id,
            'group': data.group_id,
            'database': "COG,EC,KO",
            "hsp": "mp",
            'group_detail': data.group_detail,
            "analysis_type": data.analysis_type,
            "group_method": data.group_method,
            "predict_id": str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        to_file = ["meta.export_otu_table_by_detail_without_blank(otu_table)",'meta.export_otu_seqs(otu_fasta)', "meta.export_group_table_by_detail(group)"]

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="PICRUSt2/" + main_table_name,
                            module_type=task_type,
                            params=params_json,
                            to_file=to_file)
        task_info = super(Picrust2Action, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        print(task_info)
        return json.dumps(task_info)
