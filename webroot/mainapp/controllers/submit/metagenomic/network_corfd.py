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


class NetworkCorfdAction(MetagenomicController):
    def __init__(self):
        super(NetworkCorfdAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['submit_location', 'group_id', 'group_detail', 'geneset_id', 'method', 'coefficient',
                        'coefficient_value', 'pvalue', 'task_type',"factor1", "factor2"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        fun_database = ["kegg", "cog", "vfdb", "ardb", "card", "cazy",'go','phi','qs','mvirdb','tcdb','pfam','cyps','probio']
        task_name = "metagenomic.report.network_corfd"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        factor_dict1 = json.loads(data.factor1)
        factor_dict2 = json.loads(data.factor2)
        name1 = factor_dict1["name"]
        name2 = factor_dict2["name"]
        for name in [name1,name2]:
            if not name in ["env","taxon","function"]:
                info = {"success": False, "info": "%s must be env,factor1,factor2" % name}
                return json.dumps(info)
        for factor in [factor_dict1,factor_dict2]:
            if factor["name"] == "taxon":
                for each in ["anno_id","level_id","name","top","color_level"]:
                    if not factor.has_key(each):
                        info = {"success": False, "info": "PARAMETERS MISSING: %s" % each}
                        return json.dumps(info)
                if int(factor["color_level"]) > int(factor["level_id"]):
                    info = {'success': False, 'info': '参数错误：NR的颜色显示水平必须高于等于其分类水平!', 'code': '', 'variables': ''}
                    return json.dumps(info)
                if factor["top"] == "" or int(factor["top"]) <= 1:
                    info = {"success": False, "info": "top tax must > 1"}
                    return json.dumps(info)
            elif factor["name"] == "function":
                for each in ["anno_id","level_id","name","top","database"]:
                    if not factor.has_key(each):
                        info = {"success": False, "info": "PARAMETERS MISSING: %s" % each}
                        return json.dumps(info)
                if factor["top"] == "" or int(factor["top"]) <= 1:
                    info = {"success": False, "info": "top funcion must > 1"}
                    return json.dumps(info)
            elif factor["name"] == "env":
                for each in ["env_labs","env_id"]:
                    if not factor.has_key(each):
                        info = {"success": False, "info": "PARAMETERS MISSING: %s" % each}
                        return json.dumps(info)
        task_id = geneset_info["task_id"]
        project_sn = geneset_info["project_sn"]
        [gene_profile, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        if float(data.coefficient_value) > 1 or float(data.coefficient_value) < 0:
            info = {"success": False, "info": "参数错误：相关系数阈值范围应该在[0-1]", 'code': '', 'variables': ''}
            return json.dumps(info)
        if float(data.pvalue) > 1 or float(data.pvalue) <= 0:
            info = {"success": False, "info": "参数错误：P-value范围应该在（0,1]", 'code': '', 'variables': ''}
            return json.dumps(info)
        to_file = []
        group_detail = group_detail_sort(data.group_detail)
        group_id = "all" if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "geneset_id": data.geneset_id,
            "method": data.method,
            "coefficient": data.coefficient,
            "coefficient_value": float(data.coefficient_value),
            "pvalue": float(data.pvalue),
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("geneset_id", ObjectId(data.geneset_id)),
            ("status", "start"),
            #("name", name),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            #("specimen", ",".join(samples_ids)),
        ]
        options = {
            "gene_profile": gene_profile,
            "gene_list": gene_list,
            "coefficient": data.coefficient,
            "coefficient_value": data.coefficient_value,
            "pvalue": data.pvalue,
        }
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_networkflow_personal(data,options,params_json,gene_profile,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        params_json, options = self.factor_parse(factor_dict1, factor_dict2, params_json, options)
        fac1_name =  "_" + factor_dict1["name"].upper()
        fac2_name =  "_" + factor_dict2["name"].upper()
        anno_level_name = fac1_name[0:4] + fac2_name[0:4]
        name = "CorrNetwork" + anno_level_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        name = name.replace(" ","_")
        mongo_data.append(('name', name))
        if factor_dict1["name"] == "env" or factor_dict2["name"] == "env":
            to_file.append('metagenomic.export_float_env(env_file)')
        options['group_table'] = data.geneset_id
        options["group_detail"] = data.group_detail
        to_file.append('metagenomic.export_group_table_by_detail(group_table)')
        # options['params'] = json.dumps(params_json, sort_keys=True, separators=(",", ":"))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metagenomic.insert_main_table("network_corfd", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "network_corfd"}
        options["update_info"] = json.dumps(update_info)
        #options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name="CorrNetwork/" + name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(NetworkCorfdAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

    def add_fun_args(self, fac_type, factor_dict, params_json, options):
        if fac_type == "fac1":
            params_json["factor1"] = factor_dict
        else:
            params_json["factor2"] = factor_dict
        anno_id = factor_dict["anno_id"]
        if factor_dict["name"] == "function":
            database = factor_dict["database"]
            anno_collection = "anno_" + database
            options[fac_type + "_database"] = database
        else:
            anno_collection = "anno_nr"
            options[fac_type + "_database"] = "anno_nr"
        metagenomic = Metagenomic()
        anno_info = metagenomic.from_id_get_result(anno_collection, anno_id)
        anno_file = anno_info["anno_file"]
        print anno_file
        options[fac_type + "_anno"] = anno_file
        options[fac_type + "_level"] = self.level_id(factor_dict["level_id"])
        options[fac_type + "_top"] = factor_dict["top"]
        if "color_level" in factor_dict.keys():
            options[fac_type + "_color_level"] = self.level_id(factor_dict["color_level"])
        if "second_level" in factor_dict.keys() and factor_dict["second_level"] != "":
            options[fac_type + "_second_level"] = self.level_convert(factor_dict["second_level"],factor_dict["level_id"])
            options[fac_type + "_lowestlevel"] = anno_info['lowest_level']

    def add_env_args(self, fac_type, factor_dict, params_json, options):
        if fac_type == "fac1":
            params_json["factor1"] = factor_dict
        else:
            params_json["factor2"] = factor_dict
        options['env_file'] = factor_dict["env_id"]
        options['env_id'] = factor_dict["env_id"]
        options['env_labs'] = factor_dict["env_labs"]
        return params_json, options

    def factor_parse(self,factor_dict1,factor_dict2, params_json, options):
        #factor_dict1 = json.loads(factor_dict1)
        name1 = factor_dict1["name"]
        #factor_dict2 = json.loads(factor_dict2)
        name2 = factor_dict2["name"]
        if name1 == "env":
            self.add_env_args("fac2",factor_dict1, params_json, options)
            self.add_fun_args("fac1",factor_dict2, params_json, options)
        elif name2 == "env":
            self.add_fun_args("fac1",factor_dict1, params_json, options)
            self.add_env_args("fac2",factor_dict2, params_json, options)
        else:
            self.add_fun_args("fac1", factor_dict1, params_json, options)
            self.add_fun_args("fac2", factor_dict2, params_json, options)
        return params_json, options
