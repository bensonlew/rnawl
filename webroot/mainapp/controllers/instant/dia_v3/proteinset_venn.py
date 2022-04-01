# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.dia_controller import DiaController
from mainapp.libs.signature import check_sig
import unittest
import os

class ProteinsetVennAction(DiaController):
    def __init__(self):
        super(ProteinsetVennAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        print("lxjverynice")
        print(dict(data))
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['proteinset_id']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        sg_task = self.dia.get_task_info(data.task_id)
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
            status="end",
            version='v3',
        )
        main_id = self.dia.insert_main_table('sg_proteinset_venn', main_info)
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
        task_name = 'dia_v3.report.proteinset_venn'
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
            'info': "成功存储参数，可进行venn分析"
        }
        # 更新基因集的使用信息
        self.dia.insert_proteinset_info(data.proteinset_id, name, str(main_id))
        # self.labelfree.insert_proteinset_info(data.proteinset_id, name, str(main_id))

        #
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.dia.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.dia.update_group_compare_is_use(data.task_id, data.control_id)
        ##
        
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
        cmd += "i/labelfree/proteinset_venn "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="labelfree",
            task_type="1",
            submit_location="Proteinsetvenn",
            proteinset_id="5aec21f773e5953018000029,5aec024da4e1af1097b784eb",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
