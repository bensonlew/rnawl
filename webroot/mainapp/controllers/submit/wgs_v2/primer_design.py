# -*- coding: utf-8 -*-
# __author__ = 'Liuwentian'
# modified 20190319

import os
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController
from bson.objectid import ObjectId
from biocluster.config import Config


class PrimerDesignAction(DnaController):
    """
     引物设计接口
     marker_type：标记类型snp/indel,ssr
     data_type：选择标记位置信息 vcf,location(位置文件),custom(自定义)
     marker_detail：vcf从sg_variant_site_table取结果vcf文件，
     选择标记位置信息从sg_marker_position_table中取 location:[{location:chr:start-end},{}]
     custom:[chr:start-end,....]
     condition：设置引物设计条件 {snp:{tm1:,tm2:....},indel"{}}
    """
    def __init__(self):
        super(PrimerDesignAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params1 = ["task_id", "task_type", "submit_location", "chongmingming_result", "marker_type", "data_type",
                   "marker_detail", "condition"]
        for param in params1:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        db = data.project_type = "dna_wgs_v2"
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        vcf = ""
        try:
            project_sn = result["project_sn"]
            genome_version_id = result["genome_version_id"]
            member_id = result["member_id"]
            if data.data_type == "location" or data.data_type == "custom":
                try:
                    vcf = result["pop_final_vcf"]
                except:
                    info = {"success": False, "info": "sg_task表里没有pop_final_vcf信息，请检查!"}
                    return json.dumps(info)
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn,member_id,genome_version_id信息，请检查!"}
            return json.dumps(info)
        if data.data_type == "vcf":
            result1 = Dna(db).find_one(collection="sg_variant_site_table", query_dic={"_id": ObjectId(data.marker_detail)})
            try:
                vcf = result1["vcf_path"]
            except:
                info = {"success": False, "info": "sg_variant_site_table表里没有vcf_path信息，请检查!"}
                return json.dumps(info)
        elif data.data_type == "location":
            marker_detail = []
            result2 = Dna(db).find_many(collection="sg_marker_position_table_detail",
                                       query_dic={"marker_id": ObjectId(data.marker_detail)})
            for l in result2:
                try:
                    location = l["location"]
                    marker_detail.append({'location': location})
                except:
                    info = {"success": False, "info": "sg_marker_position_table_detail中没有location信息，请检查!"}
                    return json.dumps(info)
            data.marker_detail = marker_detail
        ref_result = Dna(db).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        try:
            ref_fa = ref_result['ref']
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref信息，请检查!"}
            return json.dumps(info)
        params_json = {
            "marker_type": data.marker_type,
            "data_type": data.data_type,
            "marker_detail": json.loads(data.marker_detail) if data.data_type == "custom" else data.marker_detail,
            "condition": json.loads(data.condition),
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "task_id": data.task_id,
            "chongmingming_result": data.chongmingming_result
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = Dna(db).set_main_table_name("Primer_Design_", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Primer_Design主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_primer", data=mongo_data)
        Dna(db).update_db_record(collection="sg_primer", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_primer"}
        options = {
            "ref_fa": ref_fa,
            "marker_type": data.marker_type,
            "data_type": data.data_type,
            "condition": data.condition,
            "vcf": vcf,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if data.data_type == "custom" and data.marker_type == "ssr":
            try:
                ssr_path = ref_result['ssr_path']
                ssr = Config().SOFTWARE_DIR + "/database/dna_geneome/" + ssr_path + "ssr.ref.result.xls"
                options["ssr"] = ssr
            except:
                info = {"success": False, "info": "sg_species_version表里没有ssr_path信息，请检查!"}
                return json.dumps(info)
        if data.data_type == "custom":
            options["marker_detail"] = data.marker_detail
        elif data.data_type == "location":
            options["marker_detail"] = json.dumps(data.marker_detail)
        print options
        main_table_name = "Primer_Design_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs_v2.report.primer_design", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="workflow", params=params,
                            db_type=db)
        task_info = super(PrimerDesignAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)
