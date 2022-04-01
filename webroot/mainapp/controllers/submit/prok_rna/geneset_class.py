# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
import unittest
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import os


class GenesetClassAction(ProkRNAController):
    def __init__(self):
        super(GenesetClassAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'submit_location', 'anno_type', 'type', 'task_type']
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "lack argument: {}".format(arg)}
                return json.dumps(info)

        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0
        num = data.geneset_id.split(",")

        task_name = 'prok_rna.report.geneset_class'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "anno_type": data.anno_type,
            "task_id": data.task_id,
            "type": data.type,
        }
        # 判断传入的基因集id是否存在
        #print(data.geneset_id)
        geneset_info = {}
        geneset_name = list()
        for gd in data.geneset_id.split(","):
            geneset_info = self.prok_rna.get_main_info(gd, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset not found"}
                return json.dumps(info)
            geneset_name.append(geneset_info['name'])
        task_info = self.prok_rna.get_task_info(geneset_info['task_id'])

        if data.anno_type == "go":
            table_name = "Go"
            collection_name = "sg_geneset_go_class"
            to_file = [
                'prok_rna.export_go_class(geneset_go)'
            ]

            option = {
                "geneset_go": data.geneset_id
            }
            if task_info.get("version", None) >= "v3":
                to_file.append('prok_rna.export_go_detail(gene_go)')
                option.update({"gene_go": data.geneset_id})

        elif data.anno_type == "cog":
            table_name = "Cog"
            collection_name = "sg_geneset_cog_class"
            to_file = 'prok_rna.export_cog_class(geneset_cog)'
            option = {"geneset_cog": data.geneset_id}
        else:
            info = {'success': False, 'info': '不支持的功能分类!'}
            return json.dumps(info)

        main_table_name = 'Geneset' + table_name + "Class_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        #print(main_table_name)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ("submit_type", submit_type)
        ]
        main_table_id = self.prok_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "anno_type": data.anno_type,
            "version": task_info.get("version", None),
            "type": data.type,
            "geneset_name": ','.join(geneset_name),
            }
        options.update(option)
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])
        task_info = super(GenesetClassAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
        geneset_info = self.prok_rna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.prok_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.prok_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/prok_rna/geneset_class "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna",
            task_type="2",
            submit_location="genesetgo",
            geneset_id="5b70dcb5a4e1af118c63916c,5b70dd58a4e1af13f89b4362",
            type="origin",
            anno_type="cog"
        )
        '''
        args = dict(
            task_id="prok_rna",
            task_type="2",
            submit_location="genesetgo",
            geneset_id="5a9de1715b8c915efb2a3293,5aa214865b8c915efb2fe5b8",
            type="origin",
            anno_type="cog"
        )
        '''
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
