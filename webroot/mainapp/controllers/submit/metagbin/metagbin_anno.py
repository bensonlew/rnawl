# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.metagbin_controller import MetagbinController
from mainapp.models.mongo.metagbin import Metagbin
from mainapp.libs.signature import check_sig
from bson import SON


class MetagbinAnnoAction(MetagbinController):
    def __init__(self):
        super(MetagbinAnnoAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'genome_id', 'submit_location', 'task_type','gene_seq','gene_gff']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'metagbin.report.metagbin_anno'
        module_type = 'workflow'
        project_sn = self.metagbin.get_projectsn(data.task_id)
        update_info = {}
        params = {
                'genome_id': data.genome_id,
                'task_type': int(data.task_type),
        }
        mongo_data = [
                ('project_sn', project_sn),
                ('status', 'start'),
                ('task_id', data.task_id),
                ('desc', '正在计算'),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ]
        if hasattr(data, "database"):
            mongo_data.append(('database', data.database))
        main_table_name = 'AnnoStat_' + data.genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        params['submit_location'] = "anno_stat"
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        mongo_data.append(("name",main_table_name))
        mongo_data.append(("params", params))
        main_id = self.metagbin.insert_main_table('anno_stat', mongo_data)
        Metagbin().insert_main_table_new("anno_stat", str(main_id), {"main_id": main_id})
        update_info[str(main_id)] = 'anno_stat'
        options = {
            'sample': data.genome_id,
            'update_info': json.dumps(update_info),
            "gene_seq": data.gene_seq,
            "gene_gff": data.gene_gff,
            "main_id": str(main_id),
            }
        if hasattr(data, "database"):
            options["database"] = str(data.database)
        self.set_sheet_data(
            name=task_name,
            options=options,
            main_table_name="Annotation/" + main_table_name,
            module_type=module_type,
            project_sn=project_sn,
            task_id=data.task_id,
            )
        task_info = super(MetagbinAnnoAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)