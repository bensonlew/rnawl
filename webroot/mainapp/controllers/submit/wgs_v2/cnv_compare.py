# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20190307

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from bson.objectid import ObjectId
from mainapp.controllers.project.dna_controller import DnaController


class CnvCompareAction(DnaController):
    """
    cnv比较分析接口
    """
    def __init__(self):
        super(CnvCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["analysis_model", "region", "region_type", "sample", "marktype", "task_id", "task_type",
                  "submit_location"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        if data.analysis_model not in ["single", "multiple"]:
            info = {"success": False, "info": "分析模式%s必须是single or multiple" % data.analysis_model}
            return json.dumps(info)
        if data.region_type not in ["allregion", "location", "custom"]:
            info = {"success": False, "info": "基因组区域必须是allregion or location or custom" % data.analysis_model}
            return json.dumps(info)
        if data.analysis_model == 'multiple':
            for m in json.loads(data.marktype):
                if m not in ['same', 'diff', 'all']:
                    info = {"success": False, "info": "比较类型%s不合法!" % m}
                    return json.dumps(info)
        else:
            if data.marktype not in ['same', 'diff', 'all']:
                info = {"success": False, "info": "比较类型%s不合法!" % data.marktype}
                return json.dumps(info)
        if not hasattr(data, "project_type"):
            db = "dna_wgs_v2"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        elif data.project_type == 'dna_wgs':
            db = 'dna_wgs'
        else:
            db = data.project_types
        region_type = "real_region"
        if data.region_type == 'location':
            region_type = 'no_real_region'
            region = self.get_location(data.region, db)
        elif data.region_type == 'allregion':
            region = 'all'
        else:
            region = ','.join(json.loads(data.region))
        params_json = {
            "sample": json.loads(data.sample),
            "region": json.loads(data.region) if data.region_type != 'allregion' else data.region,
            "marktype": json.loads(data.marktype) if data.analysis_model == 'multiple' else data.marktype,
            "analysis_model": data.analysis_model,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "region_type": data.region_type,
            "task_id": data.task_id,
            "chongmingming_result": data.chongmingming_result
        }
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        # noinspection PyBroadException
        try:
            cnv_anno_path = result['cnv_anno_path']
            project_sn = result["project_sn"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!"}
            return json.dumps(info)
        cnv_anno_path_ = Dna(db).set_file_path(data.task_id, cnv_anno_path, data.client)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        # if data.analysis_model != 'single':
        #     sam1, sam2 = json.loads(data.sample)[0].split('|')
        #     main_table_name = 'CnvCompare_{}_{}_'.format(sam1, sam2) + \
        #                       datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        # else:
        #     main_table_name = 'CnvCompare_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        main_table_name = Dna(db).set_main_table_name("CnvCompare", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "cnv差异比较分析！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("subname", []),
            ('type', data.analysis_model)
        ]
        main_id = Dna(db).insert_main_table(collection="sg_cnv_compare", data=mongo_data)
        Dna(db).update_db_record(collection="sg_cnv_compare", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_cnv_compare"}
        if data.analysis_model == 'multiple':
            _marktype = ','.join(json.loads(data.marktype))
        else:
            temp = []
            for i in range(0, len(json.loads(data.sample))):
                temp.append(data.marktype)
            _marktype = ','.join(temp)
        options = {
            "infile_path": cnv_anno_path_,
            "region": region,
            "marktype": _marktype,
            "samples": ','.join(json.loads(data.sample)),  # a|b,b|c
            "region_type": region_type,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "analysis_model": data.analysis_model
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        self.set_sheet_data(name="wgs_v2.report.cnv_compare", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name=main_table_name, options=options, params=params,
                            db_type=db)
        task_info = super(CnvCompareAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def get_location(self, _id, db):
        pos_ = []
        main_id = Dna(db).find_one_record('sg_marker_position_table', {'_id': ObjectId(_id)})
        result = Dna(db).find_records('sg_marker_position_table_detail', {'marker_id': main_id})
        if result:
            for m in result:
                pos_.append(m['location'])
        return ','.join(pos_)

