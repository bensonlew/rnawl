# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig
import unittest
import os
import re
from bson.objectid import ObjectId


class ProteinsetStringPictureAction(ItraqTmtController):
    """
    蛋白质互作网络的接口
    """

    species_list = [30611, 9598, 61853, 9593, 9606, 9544, 9483, 30608, 9601, 9478, 10141, 10020, 10090, 9986, 10116,
                    43179, 37347, 9685, 9913, 9739, 9669, 9796, 132908, 59463, 9646, 9823, 9785, 9813, 9371, 9361,
                    28377, 9031, 13735, 9103, 59729, 8049, 31033, 8090, 8083, 69293, 99883, 8128, 7955, 13616, 9258,
                    9305, 9315, 7897, 7757, 7719, 51511, 6239, 7227, 4932, 15368, 4513, 4641, 4533, 4538, 4555,
                    4558, 4577, 59689, 3702, 3711, 3847, 3694, 4081, 4113, 29760, 88036, 3218, 3055, 45157]

    def __init__(self):
        super(ProteinsetStringPictureAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['proteinset_id', 'species', 'submit_location', 'useblast']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if str(data.species) == '':
            info = {"success": False, "info": "请选择参考物种".format(data.species)}
            return json.dumps(info)

        useblast = data.useblast
        if str(data.species) == '0':
            useblast = 'yes'
        '''
        if int(data.species) not in self.species_list:
            info = {"success": False, "info": "不能进行蛋白质互作分析，因为数据库中不存在该物种的蛋白质互作组数据！"}
            return json.dumps(info)
        '''
        task_name = 'itraq_and_tmt.report.proteinset_string_picture'
        task_type = 'workflow'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "proteinset_id": data.proteinset_id,
            "species": int(data.species),
            "useblast": useblast,
            "task_id": data.task_id
        }
        if hasattr(data, 'identity'):
            params_json.update({'identity': data.identity})
        proteinset_info = self.itraq_tmt.get_main_info(data.proteinset_id, 'sg_proteinset', data.task_id)
        if not proteinset_info:
            info = {"success": False, "info": "proteinset不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.itraq_tmt.get_task_info(proteinset_info['task_id'])
        stat_info = self.itraq_tmt.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id),
                                                           type="origin")
        origin_result = stat_info["result_dir"]
        if not origin_result.endswith("/"):
            origin_result = origin_result + "/"

        main_table_name = 'ProteinsetStriPic_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('proteinset_id', ObjectId(data.proteinset_id)),
            ('desc', 'string官网爬取中...'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.itraq_tmt.insert_main_table('sg_proteinset_string_picture', mongo_data)
        update_info = {str(main_table_id): "sg_proteinset_string_picture"}

        options = {
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "proteinset_list": data.proteinset_id,
            "species": data.species,
            "origin_result": origin_result,
            "useblast": useblast,
            "task_id": data.task_id
        }
        if hasattr(data, 'identity'):
            options.update({'identity': float(data.identity)})
        to_file = "itraq_tmt.export_protein_list_pfam(proteinset_list)"
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type=task_type,
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(ProteinsetStringPictureAction, self).POST()

        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        proteinset_info = self.itraq_tmt.insert_proteinset_info(data.proteinset_id, "sg_ppinetwork", str(main_table_id))

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
        cmd += "s/itraq_and_tmt/proteinset_string_picture "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_34861",
            task_type="2",
            submit_location="proteinsetstringpicture",
            proteinset_id="5d3136e517b2bf24db9e59c1",
            species= '3702',
            useblast= 'yes',
            #species = '3702',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
