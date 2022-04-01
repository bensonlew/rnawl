# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
import web
import json
import datetime
import unittest
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig


class ProteinsetPfamAction(ItraqTmtController):
    def __init__(self):
        super(ProteinsetPfamAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['proteinset_id', 'submit_location', 'type', 'task_type']
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "lack argument: {}".format(arg)}
                return json.dumps(info)
        proteinset_id = data.proteinset_id
        # if hasattr(data, 'proteinset_id2'):
        #     proteinset_id += ',' + data.proteinset_id2

        task_name = 'itraq_and_tmt.report.proteinset_pfam'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "proteinset_id": proteinset_id,
            "anno_type": data.anno_type,
            "task_id": data.task_id,
            "type": data.type,
        }
        # 判断传入的蛋白集id是否存在
        #print(data.proteinset_id)
        proteinset_info = {}
        for gd in proteinset_id.split(","):
            proteinset_info = self.itraq_tmt.get_main_info(gd, 'sg_proteinset', data.task_id)
            if not proteinset_info:
                info = {"success": False, "info": "proteinset not found"}
                return json.dumps(info)
        task_info = self.itraq_tmt.get_task_info(proteinset_info['task_id'])

        collection_name = "sg_proteinset_pfam"
        to_file = [
            'itraq_tmt.export_pfam_info(pfam_info)',
            'itraq_tmt.export_protein_list_pfam(proteinset_files)'
        ]

        main_table_name = 'Proteinset_Pfam_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.itraq_tmt.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}
        
        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "proteinset_files": proteinset_id,
            "pfam_info": proteinset_id,
            "type": data.type
            }
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])
        task_info = super(ProteinsetPfamAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
        proteinset_info = self.itraq_tmt.insert_proteinset_info(proteinset_id, collection_name, str(main_table_id))
        return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    import os
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/itraq_and_tmt/proteinset_pfam "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="tsg_34569",
            task_type="2",
            submit_location="proteinsetpfam",
            proteinset_id1="5d121dbd17b2bf0dfa0e53bd",
            proteinset_id2="5d121dbd17b2bf0dfa0e53bb",
            type="origin",
        )
        '''
        args = dict(
            task_id="itraq_tmt",
            task_type="2",
            submit_location="proteinsetgo",
            proteinset_id="5a9de1715b8c915efb2a3293,5aa214865b8c915efb2fe5b8",
            type="origin",
            anno_type="go"
        )
        '''
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
