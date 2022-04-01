# -*- coding: utf-8 -*-
import web
import json
from bson.objectid import ObjectId
import datetime
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mbio.packages.metagenomic.id_convert import name2id


class TaxonAnnoFilterAction(MetagenomicController):
    def __init__(self):
        super(TaxonAnnoFilterAction, self).__init__(instant=False)

    def POST(self):
        data = web.input()
        print(data)
        default_argu = ["type", "tax_anno_id", "group_id",
                        "submit_location", "task_type"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        filter_is = False
        for arg in ["level_filter", "sample_filter", "abu_filter"]:
            if hasattr(data, arg):
                filter_is = True
        if not filter_is:
            info = {"success": False, "info": "MISSING fitler params "}
            return json.dumps(info)

        task_name = "metagenomic.report.taxon_anno_filter"
        module_type = "workflow"
        task_type = data.task_type
        anno_info = self.metagenomic.common_find_one("anno_" + data.type,
                                                     {"_id": ObjectId(data.tax_anno_id)})

        name = data.type.capitalize() + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "type": data.type,
            "tax_anno_id": data.tax_anno_id,
            "submit_location": data.submit_location,
            "task_type": int(task_type),
            "group_id": data.group_id,
            "group_detail": json.loads(data.group_detail),
        }
        name_id = name2id(anno_info["task_id"], type="task")
        samples = []
        for sp_list in json.loads(data.group_detail).values():
            samples += sp_list
        samples = ','.join(samples)
        mongo_data = [
            ("project_sn", anno_info["project_sn"]),
            ("task_id", anno_info["task_id"]),
            ("status", "start"),
            ("specimen", samples),
            ("name", name),
            ("desc", "processing"),
            ("tax_anno_id", ObjectId(data.tax_anno_id)),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        options = {
            "samples": samples,
            "name2id": json.dumps(name_id),
            "anno_file": anno_info["anno_file"],
            "main_col": "anno_" + data.type
        }
        if hasattr(data, "level_filter"):
            params_json["level_filter"] = json.loads(data.level_filter)
            options["level_filter"] = data.level_filter
        if hasattr(data, "sample_filter"):
            params_json["sample_filter"] = json.loads(data.sample_filter)
            options["sample_filter"] = data.sample_filter
        if hasattr(data, "abu_filter"):
            params_json["abu_filter"] = json.loads(data.abu_filter)
            options["abu_filter"] = data.abu_filter
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metagenomic.insert_main_table("anno_" + data.type, mongo_data)
        update_info = {str(main_table_id): "anno_" + data.type}
        options["update_info"] = json.dumps(update_info)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=name.strip().split("_")[0] + '/' + name,
                            module_type=module_type, project_sn=anno_info["project_sn"],
                            task_id=anno_info["task_id"], params=params_json)
        task_info = super(TaxonAnnoFilterAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

