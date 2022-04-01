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


class DeleteAction(RefGenomeDbController):
    def __init__(self):
        super(DeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        input_data = web.input()
        basic_args = ["task_id"]
        # check arg
        for arg in basic_args:
            if not hasattr(input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901801', 'variables': variables}
                return json.dumps(info)
            if arg.lower() == "null":
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "%s : is null or NULL" % arg, 'code': 'C2901802',
                        'variables': variables}
                return json.dumps(info)

        task_info = self.ref_genome_db.get_task_info(input_data.client, task_id=input_data.task_id)
        if not task_info:
            info = {'success': False, 'info': '{} not exist in sg_task table'.format(input_data.task_id)}
            return json.dumps(info)
        task_status = task_info["task_status"]
        if task_status == "sanger_update":
            info = {'success': False, 'info': 'The reference genome is already on line, cannot be deleted!'}
            return json.dumps(info)
        if task_status == "tsg_delete":
            info = {'success': False, 'info': 'The reference genome is already deleted!'}
            return json.dumps(info)
        if task_status == "updating...":
            info = {'success': False, 'info': 'The reference genome is updating, can\'t be delete!'}
            return json.dumps(info)
        # delete genome database files
        result_dir = task_info['result_dir']
        if not os.path.exists(result_dir):
            pass
        else:
            shutil.rmtree(result_dir)
        # update task_info table
        self.ref_genome_db.update_task_info(input_data.client, task_id=input_data.task_id)
        # delete ref_rna_v2 sg_genome_db table
        genome_id = task_info['genome_id']
        self.ref_genome_db.delete_genome_db_table(input_data.client, genome_id)
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        desc = "delete genome database by {}".format(input_data.task_id)
        mongo_data = dict([('task_id', input_data.task_id),
                           ('submit_location', input_data.submit_location),
                           ('time', time),
                           ('desc', desc),
                           ('version', 'v2'),
                           ])
        self.ref_genome_db.insert_main_table(input_data.client, 'sg_delete', mongo_data)
        info = {"success": True, "info": "Success deleted"}
        return json.dumps(info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/ref_genome_db/delete "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="ref_genome_db_test_3",
            task_type="1",
            submit_location="delete"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
