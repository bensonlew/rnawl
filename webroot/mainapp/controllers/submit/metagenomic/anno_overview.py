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
from biocluster.file import exists


class AnnoOverviewAction(MetagenomicController):
    def __init__(self):
        super(AnnoOverviewAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ["main_id", "select", "creat_type", "geneset", "geneset_des", "geneset_list",
                        "task_type", "submit_location"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        main_id = data.main_id
        metagenomic = Metagenomic()
        if data.creat_type == "1":
            if not hasattr(data, "database_list"):
                info = {"success": False, "info": "PARAMETERS MISSING: database_list"}
                return json.dumps(info)
            else:
                main_id_info = metagenomic.from_id_get_result("anno_overview", main_id)
                task_id = main_id_info["task_id"]
        if data.creat_type == "2":
            main_id_info = metagenomic.from_id_get_result("anno_kegg", main_id)
            task_id = main_id_info["task_id"]
            kegg_params = eval(main_id_info["params"])
            group_detail = kegg_params["group_detail"]
            group_id = kegg_params["group_id"]
            group_detail = group_detail_sort(group_detail)
            group_id = "all" if group_id in ['all', 'All', 'ALL'] else ObjectId(group_id)
            all_samples_ids = self.samples_id_name(task_id, 1)
            select_samples_ids = self.ext_samples(group_detail)
            samples_ids = list(set(select_samples_ids).intersection(set(all_samples_ids)))
            samples_id_dic = id2name(task_id, type="task")
            samples_names = ",".join([samples_id_dic[each] for each in samples_ids])
        if data.creat_type == "3" or data.creat_type == "4":
            task_id = main_id
            main_id_info = metagenomic.find_origin_geneset_info(task_id)
        task_name = "metagenomic.report.anno_overview"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type
        project_sn = main_id_info["project_sn"]
        gene_origin_id1 = metagenomic.find_origin_genesetID("geneset", task_id)
        gene_origin_id = ObjectId(gene_origin_id1)
        geneset_info = metagenomic.get_geneset_info(gene_origin_id)
        gene_anno_info = metagenomic.get_anno_path(gene_origin_id, "kegg",task_id)
        gene_anno_file = gene_anno_info["anno_file"]
        gene_len_table = geneset_info['gene_list_length']
        gene_len_table = self.use_s3("/".join(gene_len_table.split("/")[0:len(gene_len_table.split("/")) - 2]) + "/")
        target_dir = gene_len_table
        new_overview = metagenomic.common_find_one('anno_overview',{'task_id': task_id})
        if new_overview.has_key("new_sum"):
            overview_table = new_overview["new_sum"]
            if exists(overview_table):  # 检查 s3 文件是否存在，不存在则表明个性化分析未完成
                target_dir = self.use_s3("/".join(overview_table.split("/")[0:len(overview_table.split("/")) - 1]) + "/")

        gene_profile = geneset_info['reads_num']
        gene_relative_profile = geneset_info['reads_num_relative']
        if data.geneset_list != "":
            geneset_ids = data.geneset_list.split(",")
            geneset_lists = []
            for eachid in geneset_ids:
                each_geneset_info = metagenomic.get_geneset_info(eachid)
                geneset_lists.append(each_geneset_info['gene_list_length'])
            geneset_files = ",".join(geneset_lists)
        name = "GENESET_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "select": json.loads(data.select),
            "creat_type": data.creat_type,
            "submit_location": data.submit_location,
            "task_type": task_type,
            "geneset": data.geneset,
            "geneset_list": data.geneset_list
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("name", data.geneset),
            ("status", "start"),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("catalog_total_length", ""),
            ("catalog_genes", ""),
            ("catalog_average_length", ""),
            ("gene_list_length", ""),
            ("download_file", " "),
            ("type", 2),
            ("geneset_des", data.geneset_des),
        ]
        options = {
            # "params": json.dumps(params_json, sort_keys=True, separators=(",", ":")),
            "select": data.select,
            "type": int(data.creat_type),
            #"gene_length_table": self.use_s3(gene_len_table),
            "gene_length_table": gene_len_table,
            "sum_anno_table" : target_dir,
            "gene_relative_profile": self.use_s3(gene_relative_profile),
            "gene_profile": self.use_s3(gene_profile),
            "task_id": task_id,
        }
        print data.select
        print gene_relative_profile
        print gene_profile
        print task_id
        if data.creat_type == "1":
            options["database_list"] = data.database_list
            params_json["database_list"] = data.database_list
        if data.creat_type == "2":
            options["samples"] = samples_names
            options["gene_kegg_anno"] = self.use_s3(gene_anno_file)
            print gene_anno_file
            params_json["group_id"] = str(group_id)
            params_json["group_detail"] = group_detail
        if data.geneset_list != "":
            options["geneset_list"] = geneset_files
        if data.creat_type == "3":
            options["select_type"] = data.select_type
            params_json["select_type"] = data.select_type
        # options['params'] = json.dumps(params_json, sort_keys=True, separators=(",", ":"))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metagenomic.insert_main_table("geneset", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "geneset"}
        options["update_info"] = json.dumps(update_info)
        # options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name=name.strip().split("_")[0] + '/' + data.geneset,
                            module_type=module_type, project_sn=project_sn,
                            task_id=task_id, params=params_json)
        task_info = super(AnnoOverviewAction, self).POST()
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
        #samples = ",".join(samples)
        return samples
