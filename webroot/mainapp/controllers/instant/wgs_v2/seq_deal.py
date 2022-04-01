# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20190312

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class SeqDealAction(DnaController):
    """
    局部组装序列获取接口--可以获取单条也可获取多条序列
    """
    def __init__(self):
        super(SeqDealAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["main_id", "sample_region_ids", "scaffold_ids", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        params_json = {
            "sample_region_ids": data.sample_region_ids,
            "scaffold_ids": data.scaffold_ids,
            "task_type": int(data.task_type),
            "main_id": data.main_id,
            "submit_location": data.submit_location,
            "task_id": data.task_id
        }
        if not hasattr(data, "project_type"):
            db = "dna_wgs_v2"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        elif data.project_type == "dna_wgs":
            db = "dna_wgs"
        else:
            db = data.project_type
        result = Dna(db).find_one(collection="sg_assembly", query_dic={"_id": ObjectId(data.main_id)})
        # noinspection PyBroadException
        try:
            project_sn = result["project_sn"]
            seq_path = result["seq_path"]
            member_id = result['member_id']
            params__ = result['params']
        except:
            info = {"success": False, "info": "sg_assembly表里没有project_sn, seq_path, params信息，请检查!"}
            return json.dumps(info)
        seq_path = Dna(db).set_file_path(data.task_id, seq_path, data.client)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = "Get_seq_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "--"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_sequence", data=mongo_data)
        Dna(db).update_db_record(collection="sg_sequence", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_sequence"}
        options = {
            "sample_region_ids": data.sample_region_ids,
            "scaffold_ids": data.scaffold_ids,
            "seq_path": seq_path.rstrip('/') + '/',
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "possname": json.dumps(self.make_possname(json.loads(params__)['poss'], json.loads(params__)['posname']))
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        sheet_data = self.set_sheet_data(name="wgs_v2.seq_deal", member_id=member_id, project_sn=project_sn,
                                         task_id=data.task_id, main_table_name="sequence/" + main_table_name,
                                         options=options, params=params, db_type=db, module_type='tool',
                                         target_output=True, analysis_name='seqdeal')
        # print sheet_data
        task_info = super(SeqDealAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def make_possname(self, poss, posname):
        possname = {}
        pos = poss.split('|')
        psname = posname.split('|')
        for i in range(0, len(pos)):
            possname[psname[i]] = pos[i]
        return possname

