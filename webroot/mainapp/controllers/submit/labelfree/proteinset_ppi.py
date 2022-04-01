# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.labelfree_controller import LabelfreeController
from mbio.api.to_file.labelfree import *
from mainapp.libs.signature import check_sig
import unittest
import os
import re


class ProteinsetPpiAction(LabelfreeController):
    """
    蛋白质互作网络的接口
    """

    species_list = [30611, 9598, 61853, 9593, 9606, 9544, 9483, 30608, 9601, 9478, 10141, 10020, 10090, 9986, 10116,
                    43179, 37347, 9685, 9913, 9739, 9669, 9796, 132908, 59463, 9646, 9823, 9785, 9813, 9371, 9361,
                    28377, 9031, 13735, 9103, 59729, 8049, 31033, 8090, 8083, 69293, 99883, 8128, 7955, 13616, 9258,
                    9305, 9315, 7897, 7757, 7719, 51511, 6239, 7227, 4932, 15368, 4513, 4641, 4533, 4538, 4555,
                    4558, 4577, 59689, 3702, 3711, 3847, 3694, 4081, 4113, 29760, 88036, 3218, 3055, 45157]

    def __init__(self):
        super(ProteinsetPpiAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['proteinset_id', 'taxon', 'species', 'submit_location', 'combine_score', 'task_type', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if str(data.species) == '':
            info = {"success": False, "info": "请选择参考物种".format(data.species)}
            return json.dumps(info)

        if str(data.species) == '0':
            info = {"success": False, "info": "您选择的物种不存在数据库中".format(data.species)}
            return json.dumps(info)
        '''
        if int(data.species) not in self.species_list:
            info = {"success": False, "info": "不能进行蛋白质互作分析，因为数据库中不存在该物种的蛋白质互作组数据！"}
            return json.dumps(info)
        '''
        task_name = 'labelfree.report.proteinset_ppi'
        task_type = 'workflow'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "taxon": data.taxon,
            "proteinset_id": data.proteinset_id,
            "species": data.species,
            "combine_score": data.combine_score,
            "task_id": data.task_id
        }
        proteinset_info = self.labelfree.get_main_info(data.proteinset_id, 'sg_proteinset', data.task_id)
        if not proteinset_info:
            info = {"success": False, "info": "proteinset不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.labelfree.get_task_info(proteinset_info['task_id'])

        if task_info.has_key('protein_fa'):
            seq = task_info['protein_fa']
        elif task_info.has_key('params'):
            work_flow_params = dict(task_info['params'])
            fa_dir = work_flow_params['protein_fasta']
            fa_dir = fa_dir.split("||")[1]
            if re.match(r'tsanger:',fa_dir):
                seq = fa_dir.replace('tsanger:','/mnt/ilustre/tsanger-data/')
            else:
                seq = fa_dir.replace('sanger:','/mnt/ilustre/data/')
        else:
            info = {"success": False, "info": "找不到序列文件"}
            return json.dumps(info)

        main_table_name = 'ProteinsetPPI_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('proteinset_id', ObjectId(data.proteinset_id)),
            ('desc', 'ppi_network分析中...'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.labelfree.insert_main_table('sg_proteinset_ppi', mongo_data)
        update_info = {str(main_table_id): "sg_proteinset_ppi"}

        options = {
            'update_info': json.dumps(update_info),
            "ppi_id": str(main_table_id),
            "proteinset_id": data.proteinset_id,
            "proteinset_list": data.proteinset_id,
            "species": data.species,
            "seq": seq,
            "combine_score": data.combine_score
        }
        to_file = "labelfree.export_protein_list_ppi(proteinset_list)"
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type=task_type,
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(ProteinsetPpiAction, self).POST()

        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        proteinset_info = self.labelfree.insert_proteinset_info(data.proteinset_id, "sg_ppinetwork", str(main_table_id))

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
        cmd += "s/labelfree/proteinset_ppi "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="labelfree",
            task_type="2",
            submit_location="proteinsetppi",
            proteinset_id="5aa214865b8c915efb2fe5b8",
            taxon = 'Plants',
            species= '4081',
            #species = '3702',
            combine_score = '300'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
