# -*- coding: utf-8 -*-

import datetime
import json
import os
import unittest
from collections import OrderedDict

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class CernaSankeyAction(WholeTranscriptomeController):
    def __init__(self):
        super(CernaSankeyAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "task_type", "submit_location", 'relation_type']
        all_args =  basic_args + ['relation_file', "relation_string", "file_id"]
        # check arg

        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903501', 'variables': variables}
                return json.dumps(info)

        task_info_dict = self.whole_transcriptome.get_task_info(data.task_id)
        project_sn = task_info_dict["project_sn"]
        task_id = data.task_id
        params = dict()
        for arg in all_args:
            if hasattr(data, arg):
                if arg == "task_type":
                    params.update({arg: int(getattr(data, arg))})
                else:
                    params.update({arg: getattr(data, arg)})

        name = "CernaSankey" + '_'
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
            desc='cerna sankey analysis main table',
            params=params,
            status="start"
        )
        main_id = self.whole_transcriptome.insert_main_table('cerna_sankey', main_info)

        options = dict()
        for arg in all_args:
            if hasattr(data, arg):
                options.update({arg: getattr(data, arg)})

        options.update({
            "main_id": str(main_id)
        })
        options.update({
            "update_info": json.dumps({str(main_id): "cerna_sankey"})
        })

        to_files = []
        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'whole_transcriptome.report.cerna_sankey'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(CernaSankeyAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }

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
        cmd += "i/whole_transcriptome/cerna_sankey "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_36658",
            task_type=1,
            submit_location="cerna_sankey",
            relation_type="file",
            relation_file="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/whole_rna/xx",
            relation_string=""
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
