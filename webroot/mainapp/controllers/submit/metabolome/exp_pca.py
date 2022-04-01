# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0622

import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON


class ExpPcaAction(MetabolomeController):
    def __init__(self):
        super(ExpPcaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'group_id', 'group_detail', 'metab_table',
                        'transform', 'confidence']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.exp_pca'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)
        group_detail = json.loads(data.group_detail)
        select_samples = self.ext_samples(group_detail)
        for eachgroup in group_detail.keys():
            sam_len = len(group_detail[eachgroup])
            if sam_len < 2:
                info = {"success": False, "info": "PCA分析输入的样本数每组不得少于2个!", 'code':'C2300601'}
                return json.dumps(info)
        name = "ExpPca_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "metab_table": data.metab_table,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "transform": data.transform,
            "confidence": float(data.confidence),
            "task_id": task_id
        }
        if project_type == "LC":
            if hasattr(data, "table_type"):
                table_type= data.table_type
                params_json['table_type'] = table_type
                # info = {"success": False, "info": "LC项目必须输入metab_tabel阴阳离子类型!", 'code':'C2300602'}
                # return json.dumps(info)
            else:
                table_type = 'mix'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, table_type)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, table_type)
        elif project_type == "GC":
            table_type = 'pos'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)



        if not data.transform in ["UV", "Par", "Ctr", ""]:
            info = {'success': False, 'info': '数据转换错误!', 'code':'C2300603'}
            return json.dumps(info)
        # if not data.table_type in ["pos", "neg","mix"]:
        #     info = {'success': False, 'info': '表格类型错误!', 'code':'C2300604'}
        #     return json.dumps(info)
        if not 0 < float(data.confidence) < 1:
            info = {'success': False, 'info': '置信度必须在0-1之间!', 'code':'C2300605'}
            return json.dumps(info)
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('table_type',table_type)
        ]
        # prepare option for workflow
        if data.transform =="":
           data.transform = "none"
        options = {
            "metab_table": metab_table_path,
            "metab_desc": metab_desc_path,
            "mul_type": "pca",
            "confidence": data.confidence,
            "data_trans": data.transform,
            "group_detail": data.group_detail,
            "group_table": data.group_id,
            "name": name
        }
        to_file = []
        to_file.append('metabolome.export_group_by_detail2(group_table)')
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = metabolome.insert_main_table("exp_pca", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "exp_pca"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = '2.SampleComp/02.ExpPCA/' + name
        else:
            m_table_name = name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(ExpPcaAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

    def ext_samples(self, group_detail):
        samples = []
        for each in group_detail.values():
            for i in each:
                if i not in samples:
                    samples.append(i)
        return samples
