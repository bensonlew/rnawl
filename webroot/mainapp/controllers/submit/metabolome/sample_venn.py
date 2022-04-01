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



class SampleVennAction(MetabolomeController):
    def __init__(self):
        super(SampleVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['submit_location', 'task_type', 'task_id', 'group_detail','group_id','metab_table','threshold']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.sample_venn'
        module_type = 'workflow'
        task_id = data.task_id
        #project_sn, project_type, member_id = metabolome.get_project_info(task_id)
        #project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        r = metabolome.conmon_find_one('sg_task',{'task_id': task_id})
        project_sn = r['project_sn']
        project_type = r['type']
        if "save_pdf" in r and int(r["save_pdf"]) == 1:
            name = "ExpVenn_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            main_table_name = '2.SampleComp/03.ExpVenn/' + name
        else:
            name = "Venn_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            main_table_name = name.strip().split("_")[0] + '/' + name

        params_json = {
            'metab_table' : data.metab_table,
            'group_id' : data.group_id,
            'group_detail' : json.loads(data.group_detail),
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'task_id': task_id,
            'threshold': data.threshold
        }

        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('status', 'start'),
            ("name", name),
            ("desc", "Running"),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('version','3.0')
        ]
        if project_type == "LC":
            if hasattr(data, "table_type"):
                table_type= data.table_type
                params_json['table_type'] = table_type
                # info = {"success": False, "info": "LC项目必须输入metab_tabel阴阳离子类型!", 'code':'C2300602'}
                # return json.dumps(info)
            else:
                table_type = 'mix'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, table_type)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, table_type)
        elif project_type == "GC":
            table_type = 'pos'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)

        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        metabolome = Metabolome()
        main_table_id = metabolome.insert_main_table('sample_venn', mongo_data)
        metabolome.insert_main_id('sample_venn', main_table_id)
        options = {
            "metab_table": metab_table_path,
            "group_detail": data.group_detail,
            'main_table_id': str(main_table_id),
            "threshold" : data.threshold,
            "metab_desc" : metab_desc_path,
            "name": name
        }
        if "save_pdf" in r and int(r["save_pdf"]) == 1:
            options["save_pdf"] = 1
        update_info = {str(main_table_id): 'sample_venn'}
        options["update_info"] = json.dumps(update_info)
        print '{}'.format(main_table_name)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name,
                            module_type=module_type, project_sn=project_sn,
                            task_id=task_id, params=params_json)
        task_info = super(SampleVennAction, self).POST()
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