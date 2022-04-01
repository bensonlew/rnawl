# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import web
import json
import datetime
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
import os
import unittest
from bson import ObjectId
from mainapp.controllers.core.basic import Basic


class DiffGenesetDisgenetAction(MedicalTranscriptomeController):
    def __init__(self):
        super(DiffGenesetDisgenetAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_args = ['geneset_id', 'submit_location', 'task_type', 'task_id',
                        'el', 'ei', 'dpi', 'dsi', 'padjust_method', 'score']
        for args in default_args:
            if not hasattr(data, args):
                variables = []
                variables.append(args)
                info = {'success': False, 'info': '%s参数缺少!' % args, 'code': 'C2901301', 'variables': variables}
                return json.dumps(info)

        task_name = 'medical_transcriptome.report.diff_geneset_disgenet'

        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "task_id": data.task_id,
            'el': data.el,
            "ei": data.ei,
            'dpi': data.dpi,
            'dsi': data.dsi,
            'score': data.score,
            "padjust_method": data.padjust_method,
        }
        # 判断传入的基因集id是否存在
        geneset_info = self.medical_transcriptome.get_main_info(data.geneset_id, 'sg_geneset', data.task_id)
        if not geneset_info:
            info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!", 'code': 'C2901302', 'variables': ''}
            return json.dumps(info)
        task_info = self.medical_transcriptome.get_task_info(geneset_info['task_id'])
        geneset_name = geneset_info['name']

        main_table_name = geneset_name + "_DisGeNET_Enrich_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        collection_name = "sg_diff_geneset_disgenet_enrich"
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ("version", "v1"),
            ("database", "DisGeNET v7.0"),
            ("desc", "DisGeNET enrichment main table"),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.medical_transcriptome.insert_main_table(collection_name, mongo_data)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'main_table_data': main_table_data,
            'entrez_list': data.task_id,
            'gene_list': data.geneset_id,
            "main_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            'el': data.el,
            "ei": data.ei,
            'dpi': data.dpi,
            'dsi': data.dsi,
            'score': data.score,
            "padjust_method": data.padjust_method,
            'geneset_name': geneset_name,
            "update_info": json.dumps({str(main_table_id): "sg_diff_geneset_disgenet_enrich"}),
        }

        # if data.padjust_method == "":
        #     params_json.update({"padjust_method": "BH"})
        #     options.update({"padjust_method": "BH"})
        # if data.el == "":
        #     params_json.update({"el": "Definitive"})
        #     options.update({"el": "Definitive"})
        # if data.score == "":
        #     params_json.update({"score": "0.5"})
        #     options.update({"score": "0.5"})

        to_file = ['medical_transcriptome_new.geneset_disgenet.export_entrez_list(entrez_list)',
                   'medical_transcriptome_new.medical_transcriptome.export_gene_list(gene_list)']
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=main_table_name,
                            module_type="workflow",
                            to_file=to_file,
                            task_id=task_info["task_id"],
                            project_sn=task_info["project_sn"],
                            new_task_id=new_task_id)
        task_info = super(DiffGenesetDisgenetAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
        self.medical_transcriptome.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_this(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/medical_transcriptome/diff_geneset_disgenet "
            cmd += "-b http://bcl.tsg.com "
            args = dict(
                submit_location="diffgenesetdisgenet_rich",
                geneset_id="5facb20d17b2bf65566b8881",
                task_id="8ju59of6lkdojsbreb482jveeb",
                task_type="2",
                el="Definitive,Limited",
                ei="0.1",
                dpi="0.1",
                dsi="0.1",
                score="0.5",
                padjust_method="BH",
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


if __name__ == '__main__':
    unittest.main()

