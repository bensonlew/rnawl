# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 2018.08.22

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class VariantCompareAction(DnaController):
    """
     变异位点比较分析接口
    """
    def __init__(self):
        super(VariantCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        print data.analysi_type
        is_same = [] # 存放is_same的数值，并转化为数值1或者2。
        genotype = []
        params = ["location", "funtype", "efftype", "task_id", "task_type", "variant_type", "analysi_type", "submit_location", "chongmingming_result"]  # paras 中的location指的是染色体的位置。
        # 这里分别有1和2来区分样本比较和样本组比较,还有样本组和样本一起比较。
        if data.efftype == "":
            efftype = None
        else:
            efftype = data.efftype
        if data.funtype == "":
            funtype = None
        else:
            funtype = data.funtype
        if data.location == "":
            location = None
        else:
            location = data.location
        if re.match(r".+,,", data.location):
            info = {"success": False, "info": "基因组区域不能只填染色体编号!", "code": "C3200201", "variables": ""}  # 前端在传染色的时候，格式是",,"
            return json.dumps(info)
        if data.analysi_type == '1': # 样本比较
            params.extend(["sample", "is_same", "genotype", "dep"]) # 样本比较添加如下参数。其中genotype是用于展示是纯合还是杂合。
            for param in params:
                if not hasattr(data, param):
                    var = []
                    var.append(param)
                    info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3200202", "variables": var}
                    return json.dumps(info)
            for i in json.loads(data.is_same):
                if i not in ["true", "false"]:
                    var = []
                    var.append(data.i)
                    info = {"success": False, "info": "相同与不同的类型%s不合法!, 必须为true或者false" %data.i,
                            "code": "C3200203", "variables": var}
                    return json.dumps(info)
                else:
                    if i == "true":
                        is_same_value = "0"
                        is_same.append(is_same_value)
                    elif i == "false":
                        is_same_value = "1"
                        is_same.append(is_same_value)
            for i in json.loads(data.genotype):
                if i not in ["homo|homo", "heter|heter", "homo|heter", "heter|homo"]:
                    var = []
                    var.append(data.i)
                    info = {"success": False, "info": "基因型杂合纯合%s不合法!, 必须为homo或者heter" % data.i,
                            "code": "C3200204", "variables": var}
                    return json.dumps(info)
                else:
                    if i == "homo|homo":
                        genotype_value = "0|0"
                        genotype.append(genotype_value)
                    elif i == "heter|heter":
                        genotype_value = "1|1"
                        genotype.append(genotype_value)
                    elif i == "homo|heter":
                        genotype_value = "0|1"
                        genotype.append(genotype_value)
                    elif i == "heter|homo":
                        genotype_value = "1|0"
                        genotype.append(genotype_value)
            db = "dna_evolution"
            result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
            try:
                pop_final_vcf = result['pop_final_vcf']
                project_sn = result["project_sn"]
                member_id = result["member_id"]
            except:
                info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!", "code": "C3200205",
                        "variables": ""}
                return json.dumps(info)
            pop_final_vcf = Dna("dna_evolution").set_file_path(data.task_id, pop_final_vcf, data.client)
            params_json = {
                "sample": json.loads(data.sample),
                "is_same": json.loads(data.is_same),
                "funtype": data.funtype,
                "efftype": data.efftype,
                "genotype": json.loads(data.genotype),
                "dep": json.loads(data.dep),
                "location": location,
                "task_type": int(data.task_type),  # 用于表示是workflow，还是module还是tool之类的。
                "variant_type": data.variant_type,
                "project_sn": project_sn,
                "submit_location":data.submit_location,
                "chongmingming_result": data.chongmingming_result,
                "analysi_type": data.analysi_type
            }
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
            # main_table_name = "variant_compare" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            main_table_name = Dna("dna_evolution").set_main_table_name("variant_compare", data.chongmingming_result)
            mongo_data = [
                ("name", main_table_name),
                ("status", "start"),
                ("project_sn", project_sn),
                ("task_id", data.task_id),
                ("params", params),
                ("desc", "变异位点比较分析！"),
                ("type", "sample"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            ]
            collection = "sg_variant_compare"
            main_id = Dna(db).insert_main_table(collection=collection, data=mongo_data)
            Dna(db).update_db_record(collection=collection, query_dict={"_id": main_id},
                                     update_dict={"main_id": main_id})
            update_info = {str(main_id): collection}
            options = {
                "sample": data.sample,
                "is_same": json.dumps(is_same),
                "vcf_file": pop_final_vcf,
                "funtype": data.funtype,    # 这里传值的时候传反了
                "efftype": data.efftype,   # 这里传值的时候传反了
                "dep": data.dep,
                "genotype": json.dumps(genotype),
                "region": location,
                "update_info": json.dumps(update_info),
                "main_id": str(main_id),
                "task_id": data.task_id,
                "analysis_type": data.analysi_type,
                "variant_type": data.variant_type,
                "project_sn": project_sn,
            }
            self.set_sheet_data(name="dna_evolution.report.variant_compare", member_id=member_id, project_sn=project_sn,
                                task_id=data.task_id, main_table_name="variant_compare/" + main_table_name,
                                options=options, params=params, db_type=db)
            task_info = super(VariantCompareAction, self).POST()
            task_info['id'] = str(main_id)
            return json.dumps(task_info)

        elif data.analysi_type == '2':
            params.extend(["group_id", "group_dict",  "maf", "ad", "miss"])  # 传进来的group_list存放所有的group信息。
            for param in params:
                if not hasattr(data, param):
                    var = []
                    var.append(param)
                    info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3200206", "variables": var}
                    return json.dumps(info)
            if data.group_id == "":
                var = []
                var.append(data.group_id)
                info = {"success": False, "info": "%s参数为空，请核查!" % data.group_id, "code": "C3200207",
                        "variables": var}
            group_dict_t = json.loads(data.group_dict)
            group_list = [] #将传入的dict格式的文件转化为string文件，以方便传输给下一步。
            print group_dict_t
            for key in group_dict_t.keys():
                group_list.append(key + ":" + ",".join(group_dict_t[key]))
            db = "dna_evolution"
            result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
            try:
                pop_final_vcf = result['pop_final_vcf']
                project_sn = result["project_sn"]
                member_id = result["member_id"]
            except:
                info = {"success": False, "info": "sg_task表里没有project_sn, member_id,pop_final_vcf信息，请检查!",
                        "code": "C3200208", "variables": ""}
                return json.dumps(info)
            pop_final_vcf = Dna("dna_evolution").set_file_path(data.task_id, pop_final_vcf, data.client)
            params_json = {
                "group_id": data.group_id,
                "group_dict": group_dict_t,
                "location": location,
                "funtype": data.funtype,
                "efftype": data.efftype,
                "maf": json.loads(data.maf),
                "ad": json.loads(data.ad),
                "miss1": json.loads(data.miss),
                "task_type": int(data.task_type),
                "variant_type": data.variant_type,
                "project_sn": project_sn,
                "submit_location": data.submit_location,
                "chongmingming_result": data.chongmingming_result,
                "analysi_type":data.analysi_type
            }

            # main_table_name = "variant_compare" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            main_table_name = Dna("dna_evolution").set_main_table_name("variant_compare", data.chongmingming_result)
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
            mongo_data = [
                ("name", main_table_name),
                ("status", "start"),
                ("project_sn", project_sn),
                ("task_id", data.task_id),
                ("params", params),
                ("type", "group"),
                ("desc", "group比较"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            ]

            main_id = Dna(db).insert_main_table(collection="sg_variant_compare", data=mongo_data)
            Dna(db).update_db_record(collection="sg_variant_compare", query_dict={"_id": main_id},
                                     update_dict={"main_id": main_id})
            update_info = {str(main_id): "sg_variant_compare"}
            options = {
                "group": json.dumps(group_list),
                "analysis_type": data.analysi_type,
                "region": location,
                "funtype": funtype,
                "efftype": efftype,
                "maf": data.maf,
                "ad": data.ad,
                "miss": data.miss,
                "update_info": json.dumps(update_info),
                "main_id": str(main_id),
                "task_id": data.task_id,
                "vcf_file": pop_final_vcf,
                "project_sn": project_sn,
            }
            self.set_sheet_data(name="dna_evolution.report.variant_compare", member_id=member_id, project_sn=project_sn,
                                task_id=data.task_id,
                                main_table_name=main_table_name, options=options, module_type="workflow", params=params,
                                db_type=db)
            task_info = super(VariantCompareAction, self).POST()
            task_info['id'] = str(main_id)
            return json.dumps(task_info)
        if data.analysi_type == '3':  # 样本比较和样本组比较
            params.extend(["sample", "is_same", "genotype", "dep", "group_id", "group_dict", "maf", "ad", "miss"])
            for param in params:
                if not hasattr(data, param):
                    var = []
                    var.append(param)
                    info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3200209", "variables": var}
                    return json.dumps(info)
            for i in json.loads(data.is_same):
                if i not in ["true", "false"]:
                    var = []
                    var.append(data.i)
                    info = {"success": False, "info": "相同与不同的类型%s不合法!, 必须为true或者false" % data.i,
                            "code": "C3200210", "variables": var}
                    return json.dumps(info)
                else:
                    if i == "true":
                        is_same_value = "0"
                        is_same.append(is_same_value)
                    elif i == "false":
                        is_same_value = "1"
                        is_same.append(is_same_value)
            for i in json.loads(data.genotype):
                if i not in ["homo|homo", "heter|heter", "homo|heter", "heter|homo"]:
                    var = []
                    var.append(data.i)
                    info = {"success": False, "info": "基因型杂合纯合%s不合法!, 必须为homo或者heter" % data.i,
                            "code": "C3200211", "variables": var}
                    return json.dumps(info)
                else:
                    if i == "homo|homo":
                        genotype_value = "0|0"
                        genotype.append(genotype_value)
                    elif i == "heter|heter":
                        genotype_value = "1|1"
                        genotype.append(genotype_value)
                    elif i == "homo|heter":
                        genotype_value = "0|1"
                        genotype.append(genotype_value)
                    elif i == "heter|homo":
                        genotype_value = "1|0"
                        genotype.append(genotype_value)
            if data.group_id == "":
                var = []
                var.append(data.group_id)
                info = {"success": False, "info": "%s参数为空，请核查!" % data.group_id, "code": "C3200212",
                        "variables": var}
            group_dict_t = json.loads(data.group_dict)
            group_list = []  # 将传入的dict格式的文件转化为string文件，以方便传输给下一步。
            for key in group_dict_t.keys():
                group_list.append(key + ":" + ",".join(group_dict_t[key]))
            db = "dna_evolution"
            result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
            if not result:
                var = []
                var.append(data.task_id)
                info = {"success": False,
                        "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id, "code": "C3200213",
                        "variables": var}
                return json.dumps(info)
            try:
                pop_final_vcf = result['pop_final_vcf']
                project_sn = result["project_sn"]
                member_id = result["member_id"]
            except:
                info = {"success": False, "info": "sg_task表里没有pop_final_vcf, project_sn, member_id信息，请检查!", "code": "C3200214",
                        "variables": ""}
                return json.dumps(info)
            pop_final_vcf = Dna("dna_evolution").set_file_path(data.task_id, pop_final_vcf, data.client)
            params_json = {
                "group_id": data.group_id,
                "group_dict": group_dict_t,
                # "location": location,
                "location": data.location,
                "funtype": data.funtype,
                "efftype": data.efftype,
                "maf": json.loads(data.maf),
                "ad": json.loads(data.ad),
                "miss1": json.loads(data.miss),
                "sample": json.loads(data.sample),
                "is_same": json.loads(data.is_same),
                "genotype": json.loads(data.genotype),
                "dep": json.loads(data.dep),
                "task_type": int(data.task_type),  # 用于表示是workflow，还是module还是tool之类的。
                "variant_type": data.variant_type,
                "project_sn": project_sn,
                "submit_location": data.submit_location,
                "chongmingming_result": data.chongmingming_result,
                "analysi_type": data.analysi_type
            }
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
            main_table_name = Dna("dna_evolution").set_main_table_name("variant_compare", data.chongmingming_result)
            mongo_data = [
                ("name", main_table_name),
                ("status", "start"),
                ("project_sn", project_sn),
                ("task_id", data.task_id),
                ("params", params),
                ("type", "all"),
                ("desc", "group和sample比较"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            ]
            main_id = Dna(db).insert_main_table(collection="sg_variant_compare", data=mongo_data)
            Dna(db).update_db_record(collection="sg_variant_compare", query_dict={"_id": main_id},
                                     update_dict={"main_id": main_id})
            update_info = {str(main_id): "sg_variant_compare"}
            options = {
                "group": json.dumps(group_list),
                "analysis_type": data.analysi_type,
                "region": location,
                "funtype": funtype,   # 这里传值的时候传反了
                "efftype": efftype,   # 这里传值的时候传反了
                "maf": data.maf,
                "ad": data.ad,
                "miss": data.miss,
                "update_info": json.dumps(update_info),
                "main_id": str(main_id),
                "task_id": data.task_id,
                "vcf_file": pop_final_vcf,
                "sample": data.sample,
                "is_same": json.dumps(is_same),
                "dep": data.dep,
                "genotype": json.dumps(genotype),
                "variant_type": data.variant_type,
                "project_sn": project_sn,
            }
            # main_table_name = "variant_compare" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            self.set_sheet_data(name="dna_evolution.report.variant_compare", member_id=member_id, project_sn=project_sn,
                                task_id=data.task_id,
                                main_table_name=main_table_name, options=options, module_type="workflow", params=params,
                                db_type=db)
            task_info = super(VariantCompareAction, self).POST()
            task_info['id'] = str(main_id)
            return json.dumps(task_info)
