# -*- coding: utf-8 -*-
#__author__ = 'haidong.gu'
import web
import re
import os
import json
import types
import datetime
from biocluster.config import Config
from bson import SON
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.metagbin import Metagbin
from mainapp.controllers.project.metagbin_controller import MetagbinController


class IdentifSpeciesAction(MetagbinController):
    def __init__(self):
        super(IdentifSpeciesAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        """
        设置接口及其参数
        :return:
        """
        data = web.input()
        print data
        default_argu = ['task_id', 'query', 'submit_location', 'task_type', 'ref']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing: %s'%argu}
                return json.dumps(info)
        project_sn = self.metagbin.get_projectsn(data.task_id)
        task_name = 'metagbin.report.taxon_identify'
        module_type = 'workflow'
        name = "ANI_" + data.query +"_"+ datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            'query': data.query,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "ref": data.ref
        }
        if data.ref == "user_defined":
            if hasattr(data, "file_dir_id"):
                if hasattr(data, "ref_path"):
                    params_json["file_dir_id"] = data.file_dir_id
                else:
                    info = {'success': False, 'info': 'parameters missing: ref_path'}
                    return json.dumps(info)
            elif hasattr(data, "accession"):
                if hasattr(data, "accession_value"):
                    params_json["accession"] = data.accession
                    params_json["accession_value"] = data.accession_value
                else:
                    info = {'success': False, 'info': 'parameters missing: accession_value'}
                    return json.dumps(info)
            elif hasattr(data, "ref_path"):
                # 只有ref_path而没有file_dir_id
                info = {'success': False, 'info': 'parameters missing: file_dir_id'}
                return json.dumps(info)
            elif hasattr(data, "accession_value"):
                info = {'success': False, 'info': 'parameters missing: accession'}
                return json.dumps(info)
            else:
                info = {'success': False, 'info': 'parameters missing user_defined info'}
                return json.dumps(info)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        update_info = {}
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_table_id = Metagbin().insert_main_table("identif_species", mongo_data)
        Metagbin().insert_main_table_new("identif_species", str(main_table_id), {"main_id": main_table_id})
        update_info[str(main_table_id)] = 'identif_species'
        options = {
            "update_info": json.dumps(update_info),
            "main_id": str(main_table_id),
            "bin_name": data.query,
            "method": "ani",
            "query": Metagbin().get_genome(data.task_id, data.query),
            "task_id": data.task_id
            # ("ref", data.ref)
        }
        if data.ref == "user_defined" and hasattr(data, "ref_path"):
            options["ref"] = data.ref_path
        elif data.ref == "user_defined" and hasattr(data, "accession_value"):
            ref_path = os.path.join(Config().SOFTWARE_DIR, "database/NCBI_bacgenome/fna", data.accession_value + ".fna")
            options["ref"] = ref_path
        else:
            ref_list = data.ref.split(",")
            ref_path_list = []
            for ref in ref_list:
                ref_path_list.append(Metagbin().get_genome(data.task_id, ref))
            options["ref_str"] = ",".join(ref_path_list)
        self.set_sheet_data(
            name=task_name,
            options=options,
            main_table_name= name.strip().split("_")[0] + '/' + name,
            module_type=module_type,
            project_sn=project_sn,
            task_id=data.task_id,
            params=params_json
        )
        task_info = super(IdentifSpeciesAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)