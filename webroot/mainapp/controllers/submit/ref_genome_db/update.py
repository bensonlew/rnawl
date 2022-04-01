# -*- coding: utf-8 -*-
# __author__ : shicaiping
# __date__ : 20200317

import codecs
import datetime
import json
import os
import unittest
import shutil

import web
from mainapp.controllers.project.ref_genome_db_controller import RefGenomeDbController
from mainapp.libs.signature import check_sig


class UpdateAction(RefGenomeDbController):
    def __init__(self):
        super(UpdateAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        input_data = web.input()
        print(input_data)
        basic_args = ["task_id", "submit_location"]
        # check arg
        for arg in basic_args:
            if not hasattr(input_data, arg):
                variables = list()
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901801', 'variables': variables}
                return json.dumps(info)
            if arg.lower() == "null":
                variables = list()
                variables.append(arg)
                info = {'success': False, 'info': "%s : is null or NULL" % arg, 'code': 'C2901802',
                        'variables': variables}
                return json.dumps(info)

        try:
            task_info = self.ref_genome_db.get_task_info(input_data.client, task_id=input_data.task_id)
        except Exception as e:
            info = {'success': False, 'info': '{}'.format(e)}
            print("info1: {}".format(info))
            return json.dumps(info)

        if not task_info:
            info = {'success': False, 'info': '{} not exist in sg_task table'.format(input_data.task_id)}
            return json.dumps(info)
        task_status = task_info["task_status"]
        if task_status == "sanger_update":
            info = {'success': False, 'info': 'The reference genome is already on line!'}
            return json.dumps(info)
        if task_status == "tsg_delete":
            info = {'success': False, 'info': 'The reference genome has been deleted!'}
            return json.dumps(info)
        if task_status == "updating...":
            info = {'success': False, 'info': 'The reference genome is updating, please don\'t repetitive submission!'}
            return json.dumps(info)
        elif task_status not in ["tsg_finish", "sanger_finish", "isanger_finish"]:
            info = {'success': False, 'info': 'task status is not ready to update!'}
            return json.dumps(info)
        else:
            self.ref_genome_db.change_task_status(input_data.client, task_id=input_data.task_id)

        try:
            project_sn = task_info["project_sn"] if task_info["project_sn"] else "None"
            genome_id = task_info["genome_id"]
            task_status = task_info["task_status"]
            task_id = task_info["task_id"]
        except Exception as e:
            info = {'success': False, 'info': '{}'.format(e)}
            print("info2: {}".format(info))
            return json.dumps(info)

        if not genome_id:
            info = {'success': False, 'info': 'Genome id not exists!'}
            return json.dumps(info)
        result_dir = task_info["result_dir"]
        # if not os.path.exists(result_dir):
        #     info = {'success': False, 'info': 'Genome file not exists!'}
        #     return json.dumps(info)

        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=input_data.submit_location,
            task_type=int(input_data.task_type),
            project_sn=project_sn,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        desc = "update_task_id_{}".format(task_id)

        # prepare main table info
        main_info = dict(
            name=desc,
            task_id=task_id,
            project_sn=project_sn,
            desc=desc,
            status="start",
            params=params,
            created_ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            version='v2'
        )
        try:
            main_id = self.ref_genome_db.insert_main_table(input_data.client, 'sg_update', main_info)
        except Exception as e:
            info = {'success': False, 'info': '{}'.format(e)}
            print("info3: {}".format(info))
            return json.dumps(info)

        # prepare option for workflow
        try:
            options = {
                "genome_id": genome_id,
                "task_status": task_status,
                "result_dir": result_dir,
                "name": desc,
                "update_info": json.dumps({str(main_id): "sg_update"})  # to update sg_status
            }

            # 把参数交给workflow运行相应的tool
            task_name = 'ref_genome_db_v2.report.update'
            self.set_sheet_data(
                                name=task_name,
                                options=options,
                                main_table_name=desc,  # 设置交互分析结果目录名
                                module_type="workflow",
                                project_sn=project_sn,
                                task_id=task_id,
                                client=input_data.client,
            )

            # 运行workflow 并传回参数
            task_info = super(UpdateAction, self).POST()
        except Exception as e:
            info = {'success': False, 'info': '{}'.format(e)}
            print("info4: {}".format(info))
            return json.dumps(info)
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': desc
            }
        }
        # task_info['group_dict'] = group_dict
        print("info5: {}".format(task_info))
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
        cmd += "s/ref_genome_db/update "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="ref_genome_db_test_10",
            task_type="2",
            submit_location="update",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
