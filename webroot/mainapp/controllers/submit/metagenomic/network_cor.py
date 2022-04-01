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
from mbio.packages.metagenomic.id_convert import id2name
from mbio.packages.metagenomic.id_convert import name2id
from .comm_creat import CommCreatAction


class NetworkCorAction(MetagenomicController):
    def __init__(self):
        super(NetworkCorAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['submit_location', 'group_id', 'group_detail', 'geneset_id', 'anno_type', 'anno_id',
                        'method', 'level', 'top', 'coefficient', 'coefficient_value','pvalue','task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        database_type = ["nr", "kegg", "cog", "vfdb", "ardb", "card", "cazy", 'go','phi','qs','mvirdb','tcdb','pfam','cyps','probio']
        if not data.anno_type in database_type:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of anno_type (%s)' % data.anno_type}
            return json.dumps(info)
        if data.anno_type == "nr":
            if hasattr(data, "color_level"):
                if int(data.color_level) > int(data.level):
                    info = {'success': False, 'info': '参数错误：NR的颜色显示水平必须高于等于其分类水平!', 'code':'C2402001', 'variables':''}
                    return json.dumps(info)
        task_name = "metagenomic.report.network_cor"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        # 从geneset中获取task和project信息 # webroot/models/mongo/metagenomic.py
        task_id = geneset_info["task_id"]
        project_sn = geneset_info["project_sn"]
        [gene_profile, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        anno_collection = "anno_" + data.anno_type
        anno_info = metagenomic.from_id_get_result(anno_collection, data.anno_id)
        anno_file = anno_info["anno_file"]
        # if not os.path.exists(anno_file):
        #     info = {'success': False, 'info': '找不到anno_file信息!', 'code':'C2402002', 'variables':''}
        #     return json.dumps(info)
        if not int(data.top) > 1:
            info = {"success": False, "info": "参数错误：总丰度前N的物种或功能必须为大于1的整数！", 'code':'C2402003', 'variables':''}
            return json.dumps(info)
        if float(data.coefficient_value) > 1 or float(data.coefficient_value) < 0:
            info = {"success": False, "info": "参数错误：相关系数阈值范围应该在[0-1]", 'code':'C2402004', 'variables':''}
            return json.dumps(info)
        if float(data.pvalue) > 1 or float(data.pvalue) <= 0:
            info = {"success": False, "info": "参数错误：P-value范围应该在（0,1]", 'code':'C2402005', 'variables':''}
            return json.dumps(info)
        group_detail = group_detail_sort(data.group_detail)
        group_id = "all" if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        #anno_level_name =  data.coefficient.capitalize() + "_" + data.anno_type.upper() +"_" + self.level_id(data.level)
        #name = "CorrNetwork" + anno_level_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
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
            "coefficient": data.coefficient,
            "coefficient_value": float(data.coefficient_value),
            "pvalue" : float(data.pvalue),
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
            #("specimen", ",".join(samples_ids)),
        ]
        level_name = self.level_id(data.level)
        to_file = []
        options = {
            "anno_file": self.use_s3(anno_file),
            "gene_profile": gene_profile,
            "gene_list": gene_list,
            "level": level_name,
            "coefficient": data.coefficient,
            "coefficient_value": data.coefficient_value,
            "pvalue" : data.pvalue,
            "top": data.top,
            "anno_type": data.anno_type,
        }
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_networkflow_personal(data,options,params_json,gene_profile,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        if hasattr(data, "color_level"):
            params_json["color_level"] = int(data.color_level)
            options["color_level"] = self.level_id(data.color_level)
        if hasattr(data, "second_level") and data.second_level != "":
            params_json["second_level"] = data.second_level
            options["second_level"] = self.level_convert(data.second_level, data.level)
            options["lowestlevel"] = anno_info['lowest_level']
        options['group_table'] = data.geneset_id
        options["group_detail"] = data.group_detail
        to_file.append('metagenomic.export_group_table_by_detail(group_table)')
        # options['params'] = json.dumps(params_json, sort_keys=True, separators=(",", ":"))
        anno_level_name = data.anno_type.upper() +"_" + self.level_id(data.level)
        name = "CorrNetwork" + anno_level_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        name = name.replace(" ","_")
        mongo_data.append(('name', name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metagenomic.insert_main_table("network_cor", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "network_cor"}
        options["update_info"] = json.dumps(update_info)
        #options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name="CorrNetwork/" + name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(NetworkCorAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

