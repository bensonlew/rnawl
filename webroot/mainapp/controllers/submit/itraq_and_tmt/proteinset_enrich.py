# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# -*- coding: utf-8 -*-
import web
import json
import datetime
import unittest
from collections import OrderedDict
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig


class ProteinsetEnrichAction(ItraqTmtController):
    def __init__(self):
        super(ProteinsetEnrichAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['proteinset_id', 'submit_location', 'task_type', 'anno_type', 'method', "type"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info':"Lack argument: %s", "variables":[argu], "code" : "C1901101"}
                return json.dumps(info)
        
        task_name = 'itraq_and_tmt.report.proteinset_enrich'
        task_type = 'workflow'

        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "anno_type": data.anno_type,
            "method": data.method,
            "type": data.type,
            "proteinset_id": data.proteinset_id,
            "task_id": data.task_id,
        }
        # 判断传入的基因集id是否存在
        proteinset_info = self.itraq_tmt.get_main_info(data.proteinset_id, 'sg_proteinset', data.task_id)
        if not proteinset_info:
            info = {"success": False, "info": "proteinset not exist", "code" : "C1901102"}
            return json.dumps(info)
        task_info = self.itraq_tmt.get_task_info(proteinset_info['task_id'])

        task_version = "V2"
        if "version" in task_info:
            task_version = task_info["version"]

        if data.anno_type == "go":
            table_name = "Go"
            collection_name = "sg_proteinset_go_enrich"
            to_file = ['itraq_tmt.export_protein_list(proteinset_list)',
                       # 'itraq_tmt.export_all_list(all_list)',
                       'itraq_tmt.export_go_list(go_list)']
            infile = {"go_list": data.proteinset_id}
        elif data.anno_type == "kegg":
            table_name = "Kegg"
            collection_name = "sg_proteinset_kegg_enrich"
            to_file = ['itraq_tmt.export_protein_list(proteinset_list)',
                       'itraq_tmt.export_multi_protein_list(proteinset_kegg)',
                       # 'itraq_tmt.export_all_list(all_list)',
                       'itraq_tmt.export_kegg_table(kegg_table)',
                       "itraq_tmt.export_add_info(add_info)"]
            infile = {"kegg_table": data.proteinset_id, "add_info": proteinset_info['task_id']}
        else:
            info = {'success': False, 'info': 'not support this enrich type', "code" : "C1901103"}
            return json.dumps(info)

        main_table_name = 'Proteinset_enrich' + table_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

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
            "anno_type": data.anno_type,
            "method": data.method,
            "proteinset_list": data.proteinset_id,
            "type": data.type,
            "task_version": task_version
            # "all_list": data.proteinset_id
            }
        options.update(infile)
        if "database_version" in  task_info:
            go_version = task_info["database_version"].get("go", "20200628").split("_")[0]
            options.update({"go_version": go_version})

        if "database_version" in  task_info:
            kegg_version = task_info["database_version"].get("kegg", "").split("_")[0]
            if kegg_version is not None and kegg_version != "":
                options.update({"kegg_version": kegg_version})
            else:
                options.update({"kegg_version": '2017'})
        if data.anno_type == "kegg":
            options.update({"proteinset_kegg": data.proteinset_id})
            options.update({"proteinset_id": data.proteinset_id})


        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type=task_type,
                            to_file=to_file, task_id=task_info["task_id"], project_sn=task_info["project_sn"])
        task_info = super(ProteinsetEnrichAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}

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
        cmd += "s/itraq_and_tmt/proteinset_enrich "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="itraq_tmt",
            task_type="2",
            submit_location="proteinset_goenrich",
            proteinset_id="5a9de1715b8c915efb2a3293",
            type="origin",
            anno_type="go",
            method='BH'
        )
        '''
        args = dict(
            task_id="itraq_tmt",
            task_type="2",
            submit_location="proteinsetgo",
            proteinset_id="5a9de1715b8c915efb2a3293",
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

