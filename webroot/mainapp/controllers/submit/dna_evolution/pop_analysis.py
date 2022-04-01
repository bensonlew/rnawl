# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180916
# controller.submit

import web
import json
import datetime
import os
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class PopAnalysisAction(DnaController):
    """
    群体进化--群体结构三个部分的接口
    data属性中肯定有chongmingming_result字段，所以在params判断的时候，就添加上，
    如果页面没有输入chongmingming_result字段的时候，默认传入的是空字符串，生信人员只需要set_main_table_name中传
    data.chongmingming_result就行了
    """
    def __init__(self):
        super(PopAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["vcf_id", "submit_location", "group_id", "group_dict", "max_missing", "min_maf", "max_maf",
                  "mindp", "task_id", "chongmingming_result"]
        if not hasattr(data, "analysis_type"):
            info = {"success": False, "info": "缺少analysis_type参数!"}
            return json.dumps(info)
        else:
            if data.analysis_type not in ['tree', "structure", 'pca']:
                info = {"success": False, "info": "analysis_type:{}参数必须是tree or "
                                                  "structure or pca!".format(data.analysis_type)}
                return json.dumps(info)
            elif data.analysis_type == 'tree':
                params.extend(["bootstrap", "tree_type"])
            elif data.analysis_type == 'structure':
                params.extend(['kmin', 'kmax'])
            else:
                pass
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            info = {"success": False,
                    "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        chr_set = task_result['chr_set']  # 染色体与sca的个数
        vcf_file = Dna("dna_evolution").find_one(collection="sg_variant_compare_filter",
                                                 query_dic={"_id": ObjectId(data.vcf_id)})
        if not vcf_file:
            info = {"success": False,
                    "info": "sg_variant_compare_filter表里没有_id: %s对应的信息，请检查!" % data.vcf_id}
            return json.dumps(info)
        else:
            # noinspection PyBroadException
            try:
                vcf_file_path = vcf_file['vcf_path']
            except:
                info = {"success": False,
                        "info": "sg_variant_compare_filter中没有vcf_path字段，请核查！"}
                return json.dumps(info)
        vcf_file_path = Dna("dna_evolution").set_file_path(data.task_id, vcf_file_path, data.client)
        params_json = {
            "vcf_id": data.vcf_id,
            "group_id": data.group_id,
            "group_dict": json.loads(data.group_dict),
            "max_missing": data.max_missing,
            "min_maf": data.min_maf,
            "max_maf": data.max_maf,
            "mindp": data.mindp,
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            'task_type': data.task_type
        }
        options = {
            "vcf_file": vcf_file_path,
            "max_missing": float(data.max_missing),
            "min_maf": float(data.min_maf),
            "max_maf": float(data.max_maf),
            "sample_list": self.set_samples(data.group_dict)
        }
        if data.mindp:
            options.update({"minDP": int(data.mindp)})
        if hasattr(data, "maxdp") and data.maxdp:
            params_json.update({"maxdp": data.maxdp})
            options.update({"maxDP": int(data.maxdp)})
        else:
            params_json.update({"maxdp": ''})
        if data.analysis_type == 'tree':
            params_json.update({'bootstrap': data.bootstrap, "tree_type": data.tree_type})
            options.update({"bs_trees": int(data.bootstrap), "tree_type": data.tree_type})
            # main_table_name = 'PopTree_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
            main_table_name = Dna("dna_evolution").set_main_table_name("PopTree", data.chongmingming_result)
            desc = 'Tree正在分析中'
            collection = "sg_pop_tree"
        elif data.analysis_type == 'structure':
            params_json.update({"kmin": data.kmin, "kmax": data.kmax})
            options.update({"chr_set": int(chr_set), "kmin": int(data.kmin), "kmax": int(data.kmax)})
            # main_table_name = 'PopStructure_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
            main_table_name = Dna("dna_evolution").set_main_table_name("PopStructure", data.chongmingming_result)
            desc = "Structure正在分析中"
            collection = "sg_pop_structure"
        else:
            # main_table_name = 'PopPca_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
            main_table_name = Dna("dna_evolution").set_main_table_name("PopPca", data.chongmingming_result)
            options.update({"chr_set": int(chr_set)})
            desc = 'PCA正在分析中'
            collection = "sg_pop_pca"
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))     # 返回json
        mongo_data = [
            ("name", main_table_name),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("status", "start"),
            ("params", params),
            ("desc", desc)
        ]
        main_id = Dna("dna_evolution").insert_main_table(collection=collection, data=mongo_data)
        Dna("dna_evolution").update_db_record(collection=collection, query_dict={"_id": main_id},
                                              update_dict={"main_id": main_id})
        update_info = {str(main_id): collection}
        options.update({
            "main_id": str(main_id),
            "update_info": json.dumps(update_info),
            'analysis_type': data.analysis_type
        })
        print options
        self.set_sheet_data(name="dna_evolution.report.pop_analysis", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="POP/" + main_table_name, options=options,
                            module_type="workflow", params=params, db_type="dna_evolution",
                            analysis_name=data.analysis_type)
        task_info = super(PopAnalysisAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def set_samples(self, group_dict):
        """
        根据data.group_dict获取到样本的列表，并以逗号分隔
        {\"all\":[\"JY102\",\"YC_bulk\"]}
        :param group_dict:
        :return:
        """
        samples = []
        group_dict_ = json.loads(group_dict)
        for key in group_dict_:
            samples.extend(group_dict_[key])
        return ','.join(samples)

