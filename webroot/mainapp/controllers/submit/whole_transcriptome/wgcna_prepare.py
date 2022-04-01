# -*- coding: utf-8 -*-

import datetime
import json
import os
import unittest
from collections import OrderedDict

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class WgcnaPrepareAction(WholeTranscriptomeController):
    def __init__(self):
        super(WgcnaPrepareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type', 'category']
        basic_args += ['group_id', 'group_dict', 'me', 'geneset_id', 'level', 'cv']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903501', 'variables': variables}
                return json.dumps(info)
        self.exp_info = self.whole_transcriptome.get_main_info_by_record("exp", task_id=data.task_id, level=data.level,
                                                                         way="tpm")
        project_sn = self.exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        exp_level = self.exp_info["level"]
        exp_type = "RSEM"

        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_id=data.group_id,
            category=data.category,
            level=data.level,
            group_dict=group_dict,
            me=data.me,
            cv=data.cv,
            geneset_id=data.geneset_id
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        # name = "Prepare" + '_' + level + '_' + quant_method + '_' + exp_type.upper() + '_'
        name = "Prepare" + '_' + data.category + '_'
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
            exp_id=str(self.exp_info['main_id']),
            desc='wgcna pre-processing analysis main table',
            level=data.level,
            category=data.category,
            params=params,
            status="start"
        )
        main_id = self.whole_transcriptome.insert_main_table('wgcna_prepare', main_info)
        if data.geneset_id.lower() not in ["all", "refall", "none"]:
            self.whole_transcriptome.insert_geneset_info(data.geneset_id, 'wgcna_prepare', str(main_id))

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": ";".join([str(self.exp_info['main_id']), data.geneset_id, data.category, data.level]),
            "group_dict": json.dumps(group_dict),
            "main_id": str(main_id),
            "me": data.me,
            "cv": data.cv,
            "level": data.level,
            "group_id": data.group_id,
            "update_info": json.dumps({str(main_id): "wgcna_prepare"})  # to update sg_status
        }
        # prepare to file
        to_files = ["whole_transcriptome.advance.export_geneset_exp_matrix(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'whole_transcriptome.report.wgcna_prepare'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
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
            _ = self.whole_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.whole_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/whole_transcriptome/wgcna_prepare "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_36088",
            task_type="2",
            submit_location="wgcna_prepare",
            group_id='5dd5152817b2bf7328f25108',
            category="mRNA",
            level='G',
            group_dict=json.dumps({"A": ["A1", "A2"], "C": ["C2", "C3"], "S": ["S1", "S3"]}).replace('"', '\\"'),
            geneset_id="5ddc80be17b2bf305f486f62",
            me="0.1",
            cv="0.1",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
