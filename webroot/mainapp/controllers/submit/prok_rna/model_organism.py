# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import web
import json
import datetime
import unittest
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import os


class ModelOrganismAction(ProkRNAController):
    def __init__(self):
        super(ModelOrganismAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['database', 'species', 'evalue', 'submit_location', 'task_id', 'task_type']

        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "parameters missing:%s" % argu}
                return json.dumps(info)

        params_json = {
            'submit_location': data.submit_location,
            'task_id': data.task_id,
            "task_type": int(data.task_type),
            'evalue': data.evalue,
            'database': data.database,
            'species': data.species,
        }

        task_info = self.prok_rna.get_task_info(data.task_id)
        if task_info['assemble_fa']:
            output_path = task_info['assemble_fa']
            cds_fa = os.path.dirname(output_path)
        else:
            info = {"success": False, "info": "该项目{}没有fasta序列文件!".format(data.task_id)}
            return info

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', 'Model Organism Annotation'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            # ('type', 'origin'),
            ('version', 'v3.1')
        ]
        if data.database == 'imodulon':
            main_table_name = 'iModulonDB'
        else:
            main_table_name = 'RegulonDB'
        main_table_name += datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data.append(('name', main_table_name))

        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.prok_rna.insert_main_table('sg_model_organism_annot', mongo_data)
        update_info = {str(main_table_id): 'sg_model_organism_annot'}

        options = {
            'update_info': json.dumps(update_info),
            'evalue': float(data.evalue),
            'cds_fa': cds_fa,
            'database': data.database,
            'species': data.species,
            'main_table_id': str(main_table_id),
            'name': main_table_name,
        }

        task_name = 'prok_rna.report.model_organism'
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name,
                            task_id=data.task_id,
                            project_sn=task_info['project_sn'],
                            module_type='workflow', params=params_json, to_file=None)

        task_info = super(ModelOrganismAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
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
        cmd += "s/prok_rna/model_organism "
        # cmd += "-b http://192.168.12.102:9090 "
        cmd += '-b http://bcl.tsg.com '
        args = dict(
            task_id="tsg_219023",
            task_type="2",
            submit_location="prokmodelorganism",
            evalue='1e-5',
            database='imodulon',
            species='e_coli',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)
if __name__ == '__main__':
    unittest.main()
