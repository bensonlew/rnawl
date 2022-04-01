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


class AnnoVfdbAction(MetagenomicController):
    def __init__(self):
        super(AnnoVfdbAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # default_argu = ["database", "geneset_id", "group_id", "group_detail", "anno_type","method", "submit_location]
        default_argu = ["database", "geneset_id", "group_id", "group_detail", "identity", "align_length", "task_type",
                        "submit_location", "align_database"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables": argu, "code" : "C2400805"}
                return json.dumps(info)
        # database_type = ["nr", "kegg", "cog", "vfdb", "ardb", "card", "cazy"]
        if not data.database == "vfdb":
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of database (%s), vfdb expected!" , "variables": data.database, "code" : "C2400806"}
            return json.dumps(info)
        print json.dumps(data, sort_keys=True, separators=(",", ":"))
        task_name = "metagenomic.report.anno_vfdb"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type
        metagenomic = Metagenomic()
        #geneset_info = metagenomic.from_id_get_result(data.geneset_id)
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        # task_info = metagenomic.get_task_info(geneset_info["task_id"])
        task_id = geneset_info["task_id"]
        project_sn = geneset_info["project_sn"]
        geneset_table = geneset_info["gene_list_length"]
        # gene_profile = geneset_info["reads_num"]
        gene_origin_id1 = metagenomic.find_origin_genesetID("geneset", task_id)
        gene_origin_id = ObjectId(gene_origin_id1)
        gene_profile = metagenomic.get_geneset_info(gene_origin_id1)["reads_num"]
        gene_anno_info = metagenomic.get_anno_path(gene_origin_id, "vfdb", task_id)
        gene_anno_file = gene_anno_info["anno_file"]
        gene_anno_dir = "/".join(gene_anno_file.split("/")[0:len(gene_anno_file.split("/")) - 1])
        if data.align_database == "Core":
            align_database = "core"
            lowest_level_file = os.path.join(gene_anno_dir, "vfdb_core_VF_profile.xls")
            gene_anno_file = os.path.join(gene_anno_dir, "gene_vfdb_core_anno.xls")
        elif data.align_database == "Predicted":
            align_database = "predict"
            lowest_level_file = os.path.join(gene_anno_dir, "vfdb_predict_VF_profile.xls")
            gene_anno_file = os.path.join(gene_anno_dir, "gene_vfdb_predict_anno.xls")
        elif data.align_database == "All":
            align_database = "all"
            lowest_level_file = os.path.join(gene_anno_dir, "vfdb_all_VF_profile.xls")
            gene_anno_file = os.path.join(gene_anno_dir, "gene_vfdb_total_anno.xls")
        group_detail = group_detail_sort(data.group_detail)
        group_id = "all" if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        select_samples_ids = self.ext_samples(group_detail)
        group_name = self.group_name(group_detail)
        all_samples_ids = self.samples_id_name(task_id, 1)
        #all_samples_names = self.samples_id_name(task_id, 2)
        samples_ids = list(set(select_samples_ids).intersection(set(all_samples_ids)))
        samples_id_dic = id2name(task_id, type="task")
        samples_names = ",".join([samples_id_dic[each] for each in samples_ids])
        name = "VFDB_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "database": data.database,
            "submit_location": data.submit_location,
            "task_type": int(task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "geneset_id": data.geneset_id,
            'align_database':data.align_database,
            # "level_filter":data.level_filter,
            "identity": 0,
            "align_length": 0
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("geneset_id", ObjectId(data.geneset_id)),
            ("status", "start"),
            ("specimen", ",".join(samples_ids)),
            ("lowest_level", "VFs"),
            ("anno_file", ""),
            ("name", name),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        options = {
            "geneset_id": data.geneset_id,
            "geneset_table": self.use_s3(geneset_table),
            "gene_anno": self.use_s3(gene_anno_file),
            "align_database": align_database,
            "lowest_level_profile": self.use_s3(lowest_level_file),
            "gene_profile": self.use_s3(gene_profile),
            "samples": samples_names,
            "database": "vfdb"
        }
        if hasattr(data, "level_filter"):
            params_json["level_filter"] = json.loads(data.level_filter)
            options["level_select"] = data.level_filter
        if hasattr(data, "sample_filter"):
            params_json["sample_filter"] = json.loads(data.sample_filter)
            sample_filter = json.loads(data.sample_filter)[0]
            mytype = sample_filter['type']
            sam_num = sample_filter['sam_num']
            sam_abu = sample_filter['abu']
            if len(samples_ids) < int(sam_num):
                info = {"success": False, "info": "样本数应少于分组方案中所选样本数!", "code":"C2400801", "variables":""}
                return json.dumps(info)
            sam_abundance = "{},{},{}".format(mytype, sam_num, sam_abu)
            options["abu_num"] = sam_abundance
        if hasattr(data, "reads_filter"):
            params_json["reads_filter"] = json.loads(data.reads_filter)
            reads_filter = json.loads(data.reads_filter)[0]
            abu_pro_type = reads_filter['type']
            abu_pro = reads_filter['proportion']
            if float(abu_pro) < 0 or float(abu_pro) > 1:
                info = {"success": False, "info": "丰度占比必须为在0-1之间的小数！", "code":"C2400802", "variables":""}
                return json.dumps(info)
            abu_proportion = str(abu_pro_type) + "," + str(abu_pro)
            options["abu_proportion"] = abu_proportion
        if float(data.identity) != 0:
            identity = float(data.identity)
            if identity < 0 or identity > 1 or len(str(data.identity).split(".")[1]) > 2:
                info = {"success": False, "info": "identity必须为在0-1之间的小数,且小数位数不超过2位！", "code":"C2400803", "variables":""}
                return json.dumps(info)
            params_json["identity"] = identity
            options['identity'] = identity
        if data.align_length != 0 or "0":
            if int(data.align_length) < 0:
                info = {"success": False, "info": "align_length必须为正整数！", "code":"C2400804", "variables":""}
                return json.dumps(info)
            params_json["align_length"] = int(data.align_length)
            options['align_length'] = int(data.align_length)
        if ObjectId(gene_origin_id1) == ObjectId(data.geneset_id):
            options['gene_type'] = "Origin"
        if group_id != "all":
            options['group'] = 2
        to_file = []
        if not group_name in ['all', 'All', 'ALL']:
            mongo_data.append(("group_name", group_name))
            options["group_detail"] = data.group_detail
            options['group_table'] = data.geneset_id
            to_file.append('metagenomic.export_group_table_by_detail(group_table)')
        else:
            mongo_data.append(("group_name", ",".join(samples_ids)))
        if hasattr(data, "funset_filter"):
            funset_filter = json.loads(data.funset_filter)
            params_json["funset_filter"] =  json.loads(data.funset_filter)
            #options['personal_funset_id'] = data.funset_filter
            options['personal_funset_file'] =  data.funset_filter
            to_file.append('metagenomic.export_function_set(personal_funset_file)')
        # options['params'] = json.dumps(params_json, sort_keys=True, separators=(",", ":"))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metagenomic.insert_main_table("anno_vfdb", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "anno_vfdb"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name=name.strip().split("_")[0] + '/' + name,
                            module_type=module_type, project_sn=project_sn,to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(AnnoVfdbAction, self).POST()
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

    def group_name(self,group_detail):
        group_name = []
        for each in group_detail.keys():
            group_name.append(each)
        return ",".join(group_name)
