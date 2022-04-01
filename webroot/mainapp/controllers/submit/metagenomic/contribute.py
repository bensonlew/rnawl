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

class ContributeAction(MetagenomicController):
    def __init__(self):
        super(ContributeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['submit_location', 'group_id', 'group_detail', 'geneset_id', 'anno_type', 'nr_anno_id',
                        "fun_anno_id", 'method', 'group_method', 'nr_level', 'fun_level', 'top_tax', 'top_fun',
                        'task_type']
        print data
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        database_type = ["nr", "kegg", "cog", "vfdb", "ardb", "card", "cazy", "gene", "annopersonal"]
        if not data.anno_type in database_type:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of anno_type (%s)' % data.anno_type}
        task_name = "metagenomic.report.contribute"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        print "geneset_info:", geneset_info
        # 从geneset中获取task和project信息 # webroot/models/mongo/metagenomic.py
        task_id = geneset_info["task_id"]
        project_sn = geneset_info["project_sn"]
        [gene_profile, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        nr_anno_info = metagenomic.from_id_get_result("anno_nr", data.nr_anno_id)
        nr_anno_file = nr_anno_info["anno_file"]
        if int(data.top_tax) <= 0 or int(data.top_tax) > 50:
            info = {"success": False, "info": "参数错误：总丰度前N的物种参数值不合法，必须为0-50！", 'code':'C2401201', 'variables':''}
            return json.dumps(info)
        if int(data.top_fun) <= 0 or int(data.top_fun) > 50:
            info = {"success": False, "info": "参数错误：总丰度前N的功能参数值不合法，必须为0-50！", 'code':'C2401202', 'variables':''}
            return json.dumps(info)
        if data.group_method == "sum":
            group_method = 1
        elif data.group_method == "average":
            group_method = 2
        elif data.group_method == "middle":
            group_method = 3
        group_detail = group_detail_sort(data.group_detail)
        group_id = "all" if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        select_samples_ids = self.ext_samples(group_detail)
        all_samples_ids = self.samples_id_name(task_id, 1)
        #all_samples_names = self.samples_id_name(task_id, 2)
        samples_ids = list(set(select_samples_ids).intersection(set(all_samples_ids)))
        samples_id_dic = id2name(task_id, type="task")
        samples_names = ",".join([samples_id_dic[each] for each in samples_ids])
        nr_level = self.level_id(data.nr_level)
        fun_level = self.level_id(data.fun_level)
        group = ",".join(group_detail.keys()) if data.group_id not in ['all', 'All', 'ALL'] else "All"
        #anno_level_name = nr_level + "_" + data.anno_type.upper() + "_" + fun_level
        #anno_level_name = data.nr_level + "_" + data.anno_type.upper() + "_" + data.fun_level
        #name = "Contribute_" + anno_level_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "anno_type": data.anno_type,
            "submit_location": data.submit_location,
            "task_type": int(task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "geneset_id": data.geneset_id,
            "method": data.method,
            "nr_anno_id": data.nr_anno_id,
            "fun_anno_id": data.fun_anno_id,
            "nr_level": int(data.nr_level),
            "fun_level": int(data.fun_level),
            "top_fun": int(data.top_fun),
            "top_tax": int(data.top_tax),
            "group_method": data.group_method
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
            ("specimen", ",".join(samples_ids)),
            ("group", group)
        ]
        options = {
            #"fun_anno_file": self.use_s3(fun_anno_file),
            "nr_anno_file": self.use_s3(nr_anno_file),
            "gene_profile": gene_profile,
            "gene_list": gene_list,
            "samples": samples_names,
            "nr_level": nr_level,
            "fun_level": fun_level,
            "top_fun": data.top_fun,
            "top_tax": data.top_tax,
            "group_method": group_method,
        }
        to_file = []
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,gene_profile,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        fun_anno_collection = "anno_" + data.anno_type
        fun_anno_info = metagenomic.from_id_get_result(fun_anno_collection, data.fun_anno_id)
        fun_anno_file = fun_anno_info["anno_file"]
        if group_id == "all":
            options["group_all"] = "T"
        options['fun_anno_file'] = self.use_s3(fun_anno_file)
        options['group_table'] = data.geneset_id
        options["group_detail"] = data.group_detail
        to_file.append('metagenomic.export_group_table_by_detail(group_table)')
        # options['params'] = json.dumps(params_json, sort_keys=True, separators=(",", ":"))
        anno_level_name = nr_level + "_" + data.anno_type.upper() + "_" + fun_level
        name = "Contribute_" + anno_level_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        name = name.replace(" ","_")
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        mongo_data.append(('name', name))
        main_table_id = self.metagenomic.insert_main_table("contribute", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "contribute"}
        options["update_info"] = json.dumps(update_info)
        #options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name=name.strip().split("_")[0] + '/' + name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(ContributeAction, self).POST()
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

    def samples_id_name(self, task_id, type):
        samples = []
        if type == 1:
            samples_dic = id2name(task_id, type="task")
        elif type == 2:
            samples_dic = name2id(task_id, type="task")
        [samples.append(i) for i in samples_dic.keys()]
        # samples = ",".join(samples)
        return samples
