# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
import os
import unittest
from bson.objectid import ObjectId


class DiffGenesetPiplineAction(MedicalTranscriptomeController):
    def __init__(self):
        super(DiffGenesetPiplineAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['diff_id', 'geneset_names']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:   %s" % arg, 'code': 'C2900301', 'variables': variables}
                return json.dumps(info)
        if str(data.geneset_names) == '[]':
            info = {'success': False, 'info': "基因集列表为空"}
            print(info)
            return json.dumps(info)
        task_info = self.medical_transcriptome.get_task_info(task_id=data.task_id)
        project_sn = task_info["project_sn"]
        task_id = data.task_id
        diff_main_table = self.db["sg_diff"].find_one({'main_id': ObjectId(data.diff_id)})
        level = diff_main_table["level"]
        species = task_info["organism_name"]
       # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            diff_id=data.diff_id,
            geneset_names=data.geneset_names,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Diff_geneset_pipline" + '_' + level  + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v1",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='differential analysis main table',
            level=level,
            geneset_names=data.geneset_names,
            params=params,
            status="start",
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_diff_geneset_pipline', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        options = {
            "diff_id": data.diff_id,
            "pipline_main_id": str(main_id),
            "level": level,
            "genesets": data.geneset_names,
            "task_type": data.task_type,
            "update_info": json.dumps({str(main_id): "sg_diff_geneset_pipline"}),
            "task_id": data.task_id,
            'species':species,
            'main_table_data': main_table_data,
        }


        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.diff_geneset_pipline'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(DiffGenesetPiplineAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
        # task_info['group_dict'] = group_dict
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
        cmd += "s/medical_transcriptome/diff_geneset_pipline "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="diffgeneset_pipline",
            diff_id='5f45c6d117b2bf78d9c9c16d',
            geneset_names=r'{"down": ["H1_vs_H3_down_12", "S1_vs_S3_down_12", "H1_vs_H2_down_12", "S1_vs_S2_down_12"], "all": ["H1_vs_H3_all_12", "S1_vs_S3_all_12", "H1_vs_H2_all_12", "S1_vs_S2_all_12"], "up": ["H1_vs_H3_up_12", "S1_vs_S3_up_12", "H1_vs_H2_up_12", "S1_vs_S2_up_12"]}'.replace('"', '\\"'),

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
