# -*- coding: utf-8 -*-
# __author__ = "qingchen.zhang"

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


class DiffKeggAction(MetagenomicController):
    def __init__(self):
        super(DiffKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # default_argu = ["database", "geneset_id", "group_id", "group_detail", "anno_type","method", "submit_location]
        default_argu = ["geneset_id", "group_id", "group_detail", "level_filter", "task_type","ci", "correction",
                        "method", "test", "tail_type", "anno_id", "submit_location"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables": argu, "code" : "C2403401"}
                return json.dumps(info)
        if float(data.ci) > 1 or float(data.ci) < 0 :
            info = {"success": False, "info": "参数错误：显著性水平参数值错误, 应该在[0,1]", "code":"C2403402", "variables":""}
            return json.dumps(info)
        if data.correction not in  ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of correction (%s)" , "variables": data.correction, "code" : "C2403403"}
            return json.dumps(info)
        if data.tail_type not in ["two.side", "greater", "less"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of tail_type (%s)" , "variables": data.tail_type, "code" : "C2403404"}
            return json.dumps(info)
        if data.test not in [ "mann", "student", "welch", "signal"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of test (%s)" , "variables": data.test, "code" : "C2403405"}
            return json.dumps(info)

        task_name = "metagenomic.report.diff_kegg"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_id = geneset_info["task_id"]
        project_sn = geneset_info["project_sn"]
        geneset_table = geneset_info["gene_list_length"]
        gene_origin_id1 = metagenomic.find_origin_genesetID("geneset", task_id)
        gene_origin_id = ObjectId(gene_origin_id1)
        #gene_profile = metagenomic.get_geneset_info(gene_origin_id1)[data.method]
        [gene_profile, gene_list] = metagenomic.export_geneset_table(gene_origin_id1, data.method)
        print "-------------%s\n-----------%s" % (gene_origin_id, task_id)
        print json.loads(data.level_filter)
        # gene_anno_info = metagenomic.get_anno_path(gene_origin_id, "kegg", task_id)
        pathway_id = json.loads(data.level_filter)[0]["name_id"]
        try:
            find_result = metagenomic.from_id_get_result("anno_kegg_pathway", data.anno_id, main_name="kegg_id", condition={"pathway": pathway_id})
        except:
            info = {"success": False, "info": "请填写正确的代谢通路ko号", "code": "C2403408"}
            return json.dumps(info)
        if len(find_result['orthology_list'].split(";")) < 3:
            info = {"success": False, "info": "too little KOs in this pathway", "code" : "C2403406"}
            return json.dumps(info)
        gene_anno_info = metagenomic.get_anno_info(data.anno_id, "kegg")
        gene_anno_file = os.path.dirname(gene_anno_info["anno_file"]) + "/"
        # xml_file = gene_anno_info["xml_file"]
        print gene_anno_file
        group_detail = group_detail_sort(data.group_detail)
        select_samples_ids,group_name_list = self.ext_samples(group_detail)
        if len(group_name_list) < 2:
            info = {"success": False, "info": "分组不可少于两组", "code" : "C2403407"}
            return json.dumps(info)
        all_samples_ids = self.samples_id_name(task_id, 1)
        samples_ids = list(set(select_samples_ids).intersection(set(all_samples_ids)))
        samples_id_dic = id2name(task_id, type="task")
        samples_names = ",".join([samples_id_dic[each] for each in samples_ids])
        name = "DIFFKEGG_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "anno_id": data.anno_id,
            "geneset_id": data.geneset_id,
            "ci": float(data.ci) if float(data.ci) != 1 else 1,
            "correction": data.correction,
            "tail_type": data.tail_type,
            "test": data.test,
            "method": data.method
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("group", ",".join(group_name_list))
        ]
        options = {
            "geneset_id": data.geneset_id,
            "group_detail": data.group_detail,
            "geneset_table": geneset_table,
            "gene_anno": gene_anno_file,
            "gene_profile": gene_profile,
            "samples": samples_names,
            "database": "kegg",
            "group_file": data.geneset_id,
            "test": data.test,
            "correction": data.correction,
            "ci": float(data.ci),
            "type": data.tail_type,
            # "xml_file": xml_file,
            "method": data.method
        }
        if hasattr(data, "level_filter"):
            params_json["level_filter"] = json.loads(data.level_filter)
            options["level_select"] = data.level_filter
        if len(set(select_samples_ids)) != len(set(all_samples_ids)):
            options['group'] = 2
        # options['params'] = json.dumps(params_json, sort_keys=True, separators=(",", ":"))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metagenomic.insert_main_table("diff_kegg", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "diff_kegg"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        if data.test == "signal":
            options['group_file'] = ':'.join([data.geneset_id, data.group_id])
            to_file = ['metagenomic.export_group_table_for_signal(group_file)']
        else:
            to_file = ['metagenomic.export_group_table_by_detail(group_file)']
        #ppm计算使用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data, options=options, params=params_json, geneset_table=gene_profile,to_file=to_file)
        #ppm计算使用
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name=name.strip().split("_")[0] + '/' + name,
                            module_type=module_type, project_sn=project_sn,
                            task_id=task_id, params=params_json, to_file=to_file)
        task_info = super(DiffKeggAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

    def ext_samples(self, group_detail):
        samples = []
        for each in group_detail.values():
            for i in each:
                if i not in samples:
                    samples.append(i)
        return samples, group_detail.keys()

    def samples_id_name(self, task_id, type):
        samples = []
        if type == 1:
            samples_dic = id2name(task_id, type="task")
        elif type == 2:
            samples_dic = name2id(task_id, type="task")
        [samples.append(i) for i in samples_dic.keys()]
        # samples = ",".join(samples)
        return samples
