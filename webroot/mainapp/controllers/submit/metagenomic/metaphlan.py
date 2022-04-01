# -*- coding: utf-8 -*
import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.models.mongo.metagenomic import Metagenomic
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mbio.packages.metagenomic.id_convert import name2id


class MetaphlanAction(MetagenomicController):
    def __init__(self):
        super(MetaphlanAction, self).__init__(instant=False)

    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['min_cu_len', "stat", "stat_q",
                        'task_id', 'in_fastq', 'qc',
                        'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        if not hasattr(data, "stat_q"):
            data.stat_q = 0.2
        if int(data.qc):
            for argu in ["qc_length", "qc_quality"]:
                if not hasattr(data, argu):
                    info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                    return json.dumps(info)
        task_name = "metagenomic.report.taxon_anno"
        module_type = "workflow"
        params_json = {
            "min_cu_len": int(data.min_cu_len),
            "stat": data.stat,
            "stat_q": float(data.stat_q),
            "task_id": data.task_id,
            "in_fastq": data.in_fastq,
            "qc": int(data.qc),
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
        }
        if int(data.qc):
            params_json.update({
                'qc_quality': int(data.qc_quality),
                "qc_length": int(data.qc_length)
            })
        task_info = self.metagenomic.get_task_info(data.task_id)
        main_table_name = "Metaphlan_Origin_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        name_id = name2id(data.task_id, type="task")
        specimen = name_id.values()
        mongo_data = {
            'project_sn': task_info['project_sn'],
            'task_id': data.task_id,
            'status': 'start',
            'desc': 'processing',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'name': main_table_name,
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'specimen': ','.join(sorted(specimen))
        }
        main_id = self.metagenomic.insert_main_table('anno_metaphlan', mongo_data)
        update_info = {str(main_id): "anno_metaphlan"}
        options = {
            "anno_meth": "metaphlan3",
            'task_id': data.task_id,
            "name2id": json.dumps(name_id),
            "min_cu_len": int(data.min_cu_len),
            "stat": data.stat,
            "stat_q": float(data.stat_q),
            "in_fastq": data.in_fastq,
            "qc": int(data.qc),
            "main_id": str(main_id),
            "main_table": "anno_metaphlan",
            "update_info": json.dumps(update_info)
        }
        if int(data.qc):
            options.update({
                'qc_quality': int(data.qc_quality),
                "qc_length": int(data.qc_length)
            })
        to_file = ["metagenomic.get_mapping_file(in_fastq)"]
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=data.task_id, project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        task_info = super(MetaphlanAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
