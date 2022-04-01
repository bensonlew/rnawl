# -*- coding: utf-8 -*-
# __author__ = 'Liuwentian'
# modified 20190328

import os
import web
import json
import datetime
from biocluster.config import Config
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController
from bson.objectid import ObjectId


class CircosAction(DnaController):
    """
     circos图接口
    """
    def __init__(self):
        super(CircosAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "submit_location", "chongmingming_result", "chromosome", "color", "variant", "color_type"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        db = data.project_type = "dna_wgs_v2"
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            genome_version_id = result["genome_version_id"]
            member_id = result["member_id"]
            pop_final_vcf = result["pop_final_vcf"]
            pop_sv_vcf = result["pop_sv_vcf"]
            cnv_anno_path = result['cnv_anno_path']
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, genome_version_id信息，请检查!"}
            return json.dumps(info)
        ref_result = Dna(db).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        try:
            ref_chrlist = ref_result['ref_chrlist']
            gff = ref_result["gff"]
            ssr_path = ref_result['ssr_path']
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref_chrlist、gff信息，请检查!"}
            return json.dumps(info)
        params_json = {
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "task_id": data.task_id,
            "chromosome": data.chromosome,
            "color": int(data.color),
            "variant": json.loads(data.variant),
            "color_type": data.color_type,
            "chongmingming_result": data.chongmingming_result
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = Dna(db).set_main_table_name("Circos", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Circos主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_circos", data=mongo_data)
        Dna(db).update_db_record(collection="sg_circos", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_circos"}
        options = {
            "gff": Config().SOFTWARE_DIR + "/database/dna_geneome/" + gff,
            "chrlist": Config().SOFTWARE_DIR + "/database/dna_geneome/" + ref_chrlist,
            "color": data.color,
            "chromosome": data.chromosome,
            # "variant": data.variant,
            "main_id": str(main_id),
            "update_info": json.dumps(update_info)
        }
        variant = json.loads(data.variant)
        # for key, value in variant.items():
        #     i = json.loads(value)
        new_variant = []
        for i in variant:
            if i["variant"] == "snp":
                if i["type"] == "before":
                    # options["snp"] = pop_final_vcf
                    i["pwd"] = pop_final_vcf
                else:
                    result = Dna(db).find_one(collection="sg_variant_site_table", query_dic={"_id": ObjectId(i["analysis_object"])})
                    try:
                        vcf_path = result["vcf_path"]
                    except:
                        info = {"success": False, "info": "sg_variant_site_table表里没有vcf_path信息，请检查!"}
                        return json.dumps(info)
                    # options["snp"] = vcf_path
                    i["pwd"] = vcf_path
            elif i["variant"] == "indel":
                if i["type"] == "before":
                    # options["indel"] = pop_final_vcf
                    i["pwd"] = pop_final_vcf
                else:
                    result = Dna(db).find_one(collection="sg_variant_site_table",
                                              query_dic={"_id": ObjectId(i["analysis_object"])})
                    try:
                        vcf_path = result["vcf_path"]
                    except:
                        info = {"success": False, "info": "sg_variant_site_table表里没有vcf_path信息，请检查!"}
                        return json.dumps(info)
                    # options["indel"] = vcf_path
                    i["pwd"] = vcf_path
            elif i["variant"] == "snpplusindel":
                if i["type"] == "before":
                    # options["snpplusindel"] = pop_final_vcf
                    i["pwd"] = pop_final_vcf
                else:
                    result = Dna(db).find_one(collection="sg_variant_site_table",
                                              query_dic={"_id": ObjectId(i["analysis_object"])})
                    try:
                        vcf_path = result["vcf_path"]
                    except:
                        info = {"success": False, "info": "sg_variant_site_table表里没有vcf_path信息，请检查!"}
                        return json.dumps(info)
                    # options["snpplusindel"] = vcf_path
                    i["pwd"] = vcf_path
            elif i["variant"] == "sv":
                if i["type"] == "after":
                    result = Dna(db).find_one(collection="sg_variant_site_table",
                                              query_dic={"_id": ObjectId(i["analysis_object"])})
                    try:
                        vcf_path = result["vcf_path"]
                    except:
                        info = {"success": False, "info": "sg_variant_site_table表里没有vcf_path信息，请检查!"}
                        return json.dumps(info)
                    # options["sv"] = vcf_path
                    i["pwd"] = vcf_path
                else:
                    # options["sv"] = pop_sv_vcf
                    i["pwd"] = pop_sv_vcf
            elif i["variant"] == "cnv":
                if i["type"] == "before":
                    cnv_anno_path_ = Dna(db).set_file_path(data.task_id, cnv_anno_path, data.client)
                    cnv_path = os.path.join(cnv_anno_path_, (i["analysis_object"] + ".cnv.anno.xls"))
                    # options["cnv"] = cnv_path
                    i["pwd"] = cnv_path
                else:
                    result = Dna(db).find_one(collection="sg_variant_site_table",
                                              query_dic={"_id": ObjectId(i["analysis_object"])})
                    try:
                        vcf_path = result["vcf_path"]
                    except:
                        info = {"success": False, "info": "sg_variant_site_table表里没有vcf_path信息，请检查!"}
                        return json.dumps(info)
                    # options["cnv"] = vcf_path
                    i["pwd"] = vcf_path
            elif i["variant"] == "ssr":
                if i["type"] == "before" or i["type"] == "ref":
                    vcf_path = Config().SOFTWARE_DIR + "/database/dna_geneome/" + ssr_path + "ssr.ref.result.xls"
                    # options["ssr"] = vcf_path
                    i["pwd"] = vcf_path
                elif i["type"] == "after":
                    result = Dna(db).find_one(collection="sg_variant_site_table",
                                              query_dic={"_id": ObjectId(i["analysis_object"])})
                    try:
                        vcf_path = result["vcf_path"]
                    except:
                        info = {"success": False, "info": "sg_variant_site_table表里没有vcf_path信息，请检查!"}
                        return json.dumps(info)
                    # options["ssr"] = vcf_path
                    i["pwd"] = vcf_path
            elif i["variant"] == "gene":
                i["pwd"] = Config().SOFTWARE_DIR + "/database/dna_geneome/" + gff
            new_variant.append(i)
            options["variant"] = json.dumps(new_variant)
        main_table_name = "Circos_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs_v2.circos", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="tool", params=params,
                            target_output=True, db_type=db)
        task_info = super(CircosAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)


# if __name__ == '__main__':
#     """
#     """
#     list = [{"variant":"gene","type":"ref","analysis_object":"","window":10000,"style":"line"},{"variant":"cnv","type":"before","analysis_object":"AH03","window":10000,"style":"heatmap"},{"variant":"indel","type":"before","analysis_object":"AH19","window":10000,"style":"scatter"},{"variant":"snp","type":"before","analysis_object":"CZ02","window":10000,"style":"line"},{"variant":"snpplusindel","type":"before","analysis_object":"AH03","window":10000,"style":"histogram"}]
#     cmd = "python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post s/wgs_v2/circos " \
#           "-c client03 -b http://192.168.12.101:9090 "
#     cmd += '-n "{}" -d "{}" '\
#         .format("task_type;chongmingming_result;submit_location;task_id;chromosome;color;variant", "1;test_liu;test_wgs;wgs_v2;chr1,chr2,chr3;2;[{\"variant\":\"gene\",\"type\":\"ref\",\"analysis_object\":\"\",\"window\":10000,\"style\":\"line\"},{\"variant\":\"cnv\",\"type\":\"before\",\"analysis_object\":\"AH03\",\"window\":10000,\"style\":\"heatmap\"},{\"variant\":\"indel\",\"type\":\"before\",\"analysis_object\":\"AH19\",\"window\":10000,\"style\":\"scatter\"},{\"variant\":\"snp\",\"type\":\"before\",\"analysis_object\":\"CZ02\",\"window\":10000,\"style\":\"line\"},{\"variant\":\"snpplusindel\",\"type\":\"before\",\"analysis_object\":\"AH03\",\"window\":10000,\"style\":\"histogram\"}]")
#     print cmd
#     os.system(cmd)