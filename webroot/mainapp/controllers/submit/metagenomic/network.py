# -*- coding: utf-8 -*-
# __author__ = "shaohua.yuan"
import web
import json
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.param_pack import group_detail_sort
import types
from mainapp.models.mongo.metagenomic import Metagenomic
import datetime
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from bson import SON
import os
from mainapp.libs.signature import check_sig
from .comm_creat import CommCreatAction

class NetworkAction(MetagenomicController):
    def __init__(self):
        super(NetworkAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['submit_location', 'group_id', 'group_detail', 'geneset_id', 'anno_type', 'anno_id',
                        'method', 'level', 'top', 'group_method', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        data.level_id = data.level # change para name for common use
        add_personal = CommCreatAction()
        info = add_personal.judge_database(data)
        if info != "":
            return json.dumps(info)
        task_name = "metagenomic.report.network"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_id = geneset_info["task_id"]
        project_sn = geneset_info["project_sn"]
        [gene_profile, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        if not int(data.top) > 0:
            info = {"success": False, "info": "参数错误：总丰度前N的物种或功能必须为正整数！", "code":"C2402101", "variables":""}
            return json.dumps(info)
        group_detail = group_detail_sort(data.group_detail)
        group_id = "all" if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        #anno_level_name = data.anno_type.upper() + "_" + self.level_id(data.level)
        #name = "Network_" + anno_level_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "anno_type": data.anno_type,
            "submit_location": data.submit_location,
            "task_type": int(task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "geneset_id": data.geneset_id,
            "method": data.method,
            "anno_id": data.anno_id,
            "level": int(data.level),
            "group_method": data.group_method,
            "top": int(data.top),
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("geneset_id", ObjectId(data.geneset_id)),
            ("status", "start"),
            #("name", name),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("anno_type", data.anno_type),
        ]
        level_name = self.level_id(data.level)
        if data.group_method == "":
            group_method = 0
        elif data.group_method == "sum":
            group_method = 1
        elif data.group_method == "average":
            group_method = 2
        elif data.group_method == "middle":
            group_method = 3
        options = {
            #"anno_file": self.use_s3(anno_file),
            "gene_profile": gene_profile,
            "gene_list": gene_list,
            "level": level_name,
            "group_method": group_method,
            "top": data.top,
            "anno_type": data.anno_type
        }
        to_file = []
        #个性化数据库测试用调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,gene_profile,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库测试用调用
        anno_collection = "anno_" + data.anno_type
        anno_info = metagenomic.from_id_get_result(anno_collection, data.anno_id)
        anno_file = anno_info["anno_file"]
        options['anno_table'] = self.use_s3(anno_file)
        options['group_table'] = data.geneset_id
        options["group_detail"] = data.group_detail
        to_file.append('metagenomic.export_group_table_by_detail(group_table)')
        # options['params'] = json.dumps(params_json, sort_keys=True, separators=(",", ":"))
        anno_level_name = data.anno_type.upper() + "_" + self.level_id(data.level)
        name = "Network_" + anno_level_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        name = name.replace(" ","_")
        mongo_data.append(('name', name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metagenomic.insert_main_table("network", mongo_data)
        update_info = {str(main_table_id): "network"}
        options["update_info"] = json.dumps(update_info)
        #options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name=name.strip().split("_")[0] + '/' + name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(NetworkAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)
