# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import web
import json
import datetime
import unittest
from mainapp.controllers.project.labelfree_controller import LabelfreeController
from mbio.api.to_file.labelfree import *
from mainapp.libs.signature import check_sig


class ProteinsetIpathAction(LabelfreeController):
    def __init__(self):
        super(ProteinsetIpathAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['proteinset_id', 'submit_location', 'type', 'task_type','task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)

        task_name = 'labelfree.report.proteinset_ipath'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "proteinset_id": data.proteinset_id,
            "task_id": data.task_id,
            "type": data.type,
        }
        # 判断传入的基因集id是否存在
        proteinset_info = {}
        for proteinset in data.proteinset_id.split(","):
            proteinset_info = self.labelfree.get_main_info(proteinset, 'sg_proteinset', data.task_id)
            if not proteinset_info:
                info = {"success": False, "info": "proteinset不存在，请确认参数是否正确！!"}
                return json.dumps(info)
        # 两个proteinset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.labelfree.get_task_info(proteinset_info['task_id'])

        table_name = "Kegg"
        collection_name = "sg_proteinset_ipath"
        to_file = ['labelfree.export_multi_protein_list(proteinset_kegg)',
                   "labelfree.export_kegg_table(kegg_table)"]

        # proteinset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个proteinset_id作为示例即可
        option = {"proteinset_kegg": data.proteinset_id,
                  "kegg_table": data.proteinset_id.split(",")[0]}

        main_table_name = 'Proteinset' + table_name + "Ipath_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            # ("pathways", ['Metabolic_pathways.svg','Regulatory_pathways.svg', 'Biosynthesis_of_secondary_metabolities.svg']),
            ("pathways", ['Metabolism.svg', 'Secondary_metabolites.svg', 'Antibiotics.svg', 'Microbial_metabolism.svg']),
        ]
        main_table_id = self.labelfree.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "proteinset_id": data.proteinset_id,
            "type":data.type,
            }
        options.update(option)
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(ProteinsetIpathAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
        proteinset_info = self.labelfree.insert_proteinset_info(data.proteinset_id, collection_name, str(main_table_id))
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
        cmd += "s/labelfree/proteinset_ipath "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="labelfree",
            task_type="2",
            submit_location="proteinsetipath",
            proteinset_id="5a9de1715b8c915efb2a3293,5aa214865b8c915efb2fe5b8",
            type="origin"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
