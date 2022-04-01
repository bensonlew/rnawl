# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class DiffGenesetDoClassAction(MedicalTranscriptomeController):
    def __init__(self):
        super(DiffGenesetDoClassAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id',  'submit_location', 'task_id', 'task_type']
        for arg in default_argu:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "lack argument: %s" % arg, 'code': 'C2901004', 'variables': variables}
                return json.dumps(info)

        task_name = 'medical_transcriptome.report.diff_geneset_do'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "task_id": data.task_id,
        }
        # 判断传入的基因集id是否存在
        # print(data.geneset_id)
        geneset_info = {}
        gn = list()
        for gd in data.geneset_id.split(","):
            geneset_info = self.medical_transcriptome.get_main_info(gd, 'sg_geneset', data.task_id)
            gn.append(geneset_info.get("name", "unknown"))
            if not geneset_info:
                info = {"success": False, "info": "geneset not found", 'code': 'C2901005', 'variables': ''}
                return json.dumps(info)
        task_info = self.medical_transcriptome.get_task_info(geneset_info['task_id'])

        source = 'diff_exp'

        table_name = "DO_class"
        collection_name = "sg_diff_geneset_do_class"
        to_file = 'medical_transcriptome_new.geneset_annot.export_do_anno_genesets(geneset_do)'

        main_table_name = table_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        # print(main_table_name)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('version', 'v1'),
            ('status', 'start'),
            ('level', 'G'),
            ('name', main_table_name),
            ('geneset_id', data.geneset_id),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.medical_transcriptome.insert_main_table(collection_name, mongo_data)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'main_table_data': main_table_data,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_do": data.geneset_id,
            "geneset_names": ",".join(gn),
            "source": source
        }
        options.update(options)
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'],
                            new_task_id=new_task_id)

        task_info = super(DiffGenesetDoClassAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}

        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
        geneset_info = self.medical_transcriptome.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))

        return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/medical_transcriptome/diff_geneset_do_class "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            submit_location="diff_geneset_do_class",
            task_type="2",
            # db_type="medical_transcriptome",
            geneset_id="5fbbdb3117b2bf2e14b048bd,5fbbdb3c17b2bf2e14b048bf",
            geneset_type="G",
            task_id="r33k8qda64pugg53k2h51mr875",

        )

        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
