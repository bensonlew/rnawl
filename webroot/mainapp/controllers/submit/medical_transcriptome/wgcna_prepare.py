# -*- coding: utf-8 -*-

import datetime
import json
import os
import unittest
from collections import OrderedDict
from bson.objectid import ObjectId
import web
from mbio.api.to_file.medical_transcriptome import *
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig


class WgcnaPrepareAction(MedicalTranscriptomeController):
    def __init__(self):
        super(WgcnaPrepareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['exp_id', 'group_id', 'group_dict', 'me', 'geneset_id', 'level', 'cv']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903501', 'variables': variables}
                return json.dumps(info)
        self.exp_info = self.medical_transcriptome.get_main_info(data.exp_id, 'sg_exp', data.task_id)
        project_sn = self.exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        level = self.exp_info["level"]
        quant_method = self.exp_info['method']
        exp_type = "RSEM"

        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_id=data.group_id,
            level=data.level,
            group_dict=group_dict,
            me=data.me,
            cv=data.cv,
            geneset_id=data.geneset_id,
            exp_id=data.exp_id
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Prepare" + '_' + level + '_' + quant_method + '_' + exp_type.upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=data.exp_id,
            desc='wgcna pre-processing analysis main table',
            level=data.level,
            params=params,
            status="start"
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_wgcna_prepare', main_info)
        if data.geneset_id.lower() not in ["all", "refall", "none"]:
            self.medical_transcriptome.insert_geneset_info(data.geneset_id, 'sg_wgcna_prepare', str(main_id))

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        connect_exp = self.db['sg_exp']
        record_exp = connect_exp.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        is_rmbe = str(record_exp['is_rmbe']).lower()
        if is_rmbe == 'false':
            exp_id = data.exp_id
        if is_rmbe == 'true':
            exp_id == str(record_exp['batch_main_id'])
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}
        options = {
            "exp_matrix": str(self.exp_info['main_id']) + ';' + data.geneset_id + ';' +
                          data.level + ';' + is_rmbe,
            "group_dict": json.dumps(group_dict),
            "main_id": str(main_id),
            "me": data.me,
            "cv": data.cv,
            "level": data.level,
            "group_id": data.group_id,
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_wgcna_prepare"})  # to update sg_status
        }
        # prepare to file
        to_files = ["medical_transcriptome.export_geneset_exp_matrix_new(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.wgcna_prepare'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(WgcnaPrepareAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }
        # task_info['group_dict'] = group_dict
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/medical_transcriptome/wgcna_prepare "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="wgcna_prepare",
            group_id='5f46228c17b2bf20e4e269e1',
            level='G',
            group_dict=json.dumps({"H1": ["H1581_1", "H1581_2", 'H1581_3'], "H2": ["H1581_4", "H1581_5", 'H1581_6'], "H3": ["H1581_7", "H1581_8", 'H1581_9'], 'S1':["SNU16_1", "SNU16_2", "SNU16_3"], 'S2':["SNU16_4", "SNU16_5", "SNU16_6"], 'S3':["SNU16_7", "SNU16_8", "SNU16_9"]}).replace('"', '\\"'),
            geneset_id="All",
            me="1",
            cv="0.1",
            exp_id='5f50cacf17b2bf5a6c8bfd88'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
