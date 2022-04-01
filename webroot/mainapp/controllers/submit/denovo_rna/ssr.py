# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
# import json
# import random
from mainapp.libs.signature import check_sig
# from bson.objectid import ObjectId
from mainapp.libs.param_pack import *
# from biocluster.config import Config
from mainapp.models.mongo.submit.denovo_rna.gene_structure import GeneStructure
# import types
from mainapp.models.mongo.meta import Meta
# from mainapp.models.workflow import Workflow
from mainapp.controllers.project.denovo_controller import DenovoController


class Ssr(DenovoController):
    def __init__(self):
        super(Ssr, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if not hasattr(data, "sequence_id"):
            info = {"success": False, "info": "缺少参数sequence_id!"}
            return json.dumps(info)
        if not hasattr(data, "orf_id"):
            info = {"success": False, "info": "缺少参数orf_id!"}
            return json.dumps(info)
        if not hasattr(data, "primer"):
            info = {"success": False, "info": "缺少参数primer!"}
            return json.dumps(info)
        #print(data.primer)
        # if data.primer == "0":
        #     primer = "true"
        # else:
        #     primer = "false"
        #print(data.sequence_id)
        main_table_name = "Ssr_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        params = {
            "sequence_id": data.sequence_id,
            "primer": data.primer,
            "orf_id": data.orf_id,
            "task_type": data.task_type,
            "submit_location": data.submit_location
        }
        # params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        sequence_info = Meta(db=self.mongodb).get_main_info(data.sequence_id, 'sg_denovo_sequence')
        if sequence_info:
            task_id = sequence_info["task_id"]
            project_sn = sequence_info["project_sn"]
            ssr_id = GeneStructure().add_ssr_table(project_sn, task_id, params=params, name=main_table_name)
            update_info = {str(ssr_id): "sg_denovo_ssr"}
            update_info = json.dumps(update_info)

            to_file = ["denovo.export_fasta_path(gene_fasta)", "denovo.export_bed_path(bed_file)"]
            options = {
                    "gene_fasta": data.sequence_id,
                    "insert_id": str(ssr_id),
                    "bed_file": data.orf_id,
                    "update_info": update_info,
                    "primer": data.primer
                }
            self.set_sheet_data(name='denovo_rna.report.ssr', options=options, main_table_name=main_table_name, module_type='workflow', to_file=to_file, main_id=ssr_id, collection_name="sg_denovo_ssr")

            task_info = super(Ssr, self).POST()
            task_info['content'] = {'ids': {'id': str(ssr_id), 'name': main_table_name}}
            #print task_info
            return json.dumps(task_info)
        else:
            info = {"success": False, "info": "sequence_id不存在，请确认参数是否正确！!"}
            return json.dumps(info)
