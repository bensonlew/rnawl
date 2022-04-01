# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mainapp.libs.signature import check_sig
import unittest
import os

class ProteinsetVennAction(ItraqTmtController):
    def __init__(self):
        super(ProteinsetVennAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['proteinset_id']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: %s", "variables":[arg], "code" : "C1901501"}
                return json.dumps(info)
        sg_task = self.itraq_tmt.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        task_id = data.task_id
        # create main table record
        # proteinset_id_list = str(data.proteinset_id).split(',')
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            proteinset_id=data.proteinset_id,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "PsetVenn"  + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Proteinset venn analysis main table',
            params=params,
            status="end"
        )
        main_id = self.itraq_tmt.insert_main_table('sg_proteinset_venn', main_info)
        collection_name = "sg_proteinset_venn"
        update_info = {str(main_id): collection_name}
        options = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_id),
            "proteinset_id": data.proteinset_id,
        }

        # 运行workflow 并传回参数
        main_table_name = 'Proteinset' + "Venn_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        task_name = 'itraq_and_tmt.report.proteinset_venn'
        to_file = ''
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=main_table_name,
                            module_type="workflow",
                            to_file=to_file,
                            project_sn=sg_task['project_sn'],
                            task_id=sg_task['task_id'])
        task_info = super(ProteinsetVennAction, self).POST()

        task_info = dict()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                },
            'info': "success"
        }
        # 更新基因集的使用信息
        self.itraq_tmt.insert_proteinset_info(data.proteinset_id, name, str(main_id))
        self.itraq_tmt.insert_proteinset_info(data.proteinset_id, name, str(main_id))

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
        cmd += "i/itraq_and_tmt/proteinset_venn "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="itraq_tmt",
            task_type="1",
            submit_location="Proteinsetvenn",
            proteinset_id="5a35da8ba4e1af38e6fb83d7,5a35da8ba4e1af38e6fb83d9,5a35da8ba4e1af38e6fb83db",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
