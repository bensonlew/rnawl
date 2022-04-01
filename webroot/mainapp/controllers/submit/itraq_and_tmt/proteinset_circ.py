# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import web
import json
import datetime
import unittest
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig


class ProteinsetCircAction(ItraqTmtController):
    def __init__(self):
        super(ProteinsetCircAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['enrich_id', 'enrich_type', 'diff_id', 'compare_group', 'submit_location', 'task_type','task_id']
        print("input par is ".format(data))
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info':"Lack argument: %s", "variables":[argu], "code" : "C1900701"}
                print("return info is {}".format(json.dumps(info)) )
                return json.dumps(info)

        task_name = 'itraq_and_tmt.report.proteinset_circ'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "enrich_id": data.enrich_id,
            "diff_id": data.diff_id,
            "compare_group": data.compare_group,
            "task_id" : data.task_id,
            "enrich_type": data.enrich_type
        }

        if data.enrich_type == "GO":
            if data.go_type not in ['All', 'ALL', 'BP', 'CC', 'MF']:
                info = {"success": False, "info": u"%s result not exist", "variables":[data.go_type], "code" : "C1900702"}
                return info
            else:
                params_json.update(dict(go_type=data.go_type.upper()))
        # 判断传入的富集id是否存在
        if data.enrich_type == "GO":
            enrich_info = self.itraq_tmt.get_main_info(data.enrich_id, 'sg_proteinset_go_enrich', data.task_id)
        elif data.enrich_type == "KEGG":
            enrich_info = self.itraq_tmt.get_main_info(data.enrich_id, 'sg_proteinset_kegg_enrich', data.task_id)
        else:
            info = {"success": False, "info": "%s result not exist", "variables":[data.enrich_type], "code" : "C1900703"}
            return info


        if not enrich_info:
            info = {"success": False, "info": "result not exist！!", "code" : "C1900704"}
            return json.dumps(info)

        task_info = self.itraq_tmt.get_task_info(enrich_info['task_id'])

        # 判断传入的差异分析是否存在 ，暂时没写

        table_name = "Kegg"
        collection_name = "sg_proteinset_circ"
        if data.enrich_type == "GO":
            to_file = ['itraq_tmt.export_go_enrich_matrix(enrich_table)',
                       "itraq_tmt.export_compare_exp_fc(diff_fc)"]
        elif data.enrich_type == "KEGG":
            to_file = ['itraq_tmt.export_kegg_enrich_matrix(enrich_table)',
                       "itraq_tmt.export_compare_exp_fc(diff_fc)"]

        # proteinset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个proteinset_id作为示例即可
        main_table_name = "Proteinset_chord_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

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
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "enrich_id": data.enrich_id,
            "diff_id": data.diff_id,
            "compare_group": data.compare_group,
            "task_id" : data.task_id,
            "enrich_type": data.enrich_type
        }

        option = {
            "diff_fc": "diff_fc",
            "enrich_table": "enrich_table"
        }
        options.update(option)

        if data.enrich_type == "GO":
            options.update(dict(go_type=data.go_type))
        else:
            pass
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(ProteinsetCircAction, self).POST()

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
        cmd += "s/itraq_and_tmt/proteinset_circ "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="itraq_tmt",
            task_type="2",
            submit_location="proteinsetcirc",
            enrich_id="5ac35a43a4e1af07a309063b",
            enrich_type="KEGG",
            compare_group="BR_3|BR_5",
            diff_id="5aa8a6a7a4e1af7c051fe187"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
