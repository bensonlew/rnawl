# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import os
import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
import sqlite3
import unittest
# __author__ = 'liulinmeng'


class MetabsetVennAction(MetabolomeController):
    def __init__(self):
        super(MetabsetVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['submit_location', 'task_type', 'task_id', 'metabset']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.metabset_venn'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, member_id = metabolome.get_project_info(task_id)
        name = "MetabsetVenn_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        params_json = {
            'metabset': data.metabset,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'task_id': task_id
        }

        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('status', 'start'),
            ("name", name),
            ("desc", "Running"),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        metabolome = Metabolome()
        main_table_id = metabolome.insert_main_table('metabset_venn', mongo_data)
        metabolome.insert_main_id('metabset_venn', main_table_id)
        options = {
            'metabset': data.metabset,
            'main_table_id': str(main_table_id),
            "name": name
        }
        update_info = {str(main_table_id): 'metabset_venn'}
        options["update_info"] = json.dumps(update_info)
        to_file = []
        to_file.append('metabolome.export_mul_metab_set(metabset)')
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = "5.Metabset/01.MetabsetVenn/" + name
        else:
            m_table_name = "Metabset/"+name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(MetabsetVennAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)



class TestFunction(unittest.TestCase):
    """
    This is test for the web_api func. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/metabolome/metabset_venn "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="metabolome",
            task_type="1",  # maybe you need to change it
            submit_location="metabset_venn",
            metabset="5b3adf5dd887f3841c000029,5b3ae73ad887f3741c000029,5b35f332d887f3781900002c",
            project_sn="project"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()