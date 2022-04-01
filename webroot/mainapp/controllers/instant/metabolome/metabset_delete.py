# -*- coding: utf-8 -*-
# __author__ = 'shijin'
import web
import json
import os
import unittest
import datetime
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.metabolome_controller import MetabolomeController


class MetabsetDeleteAction(MetabolomeController):
    def __init__(self):
        super(MetabsetDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        #print data
        default_argu = ['set_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing: %s' % argu}
                return json.dumps(info)
        #print(data.geneset_id)
        if self.Metabolome.delete_set(data.set_id):
            task_info = {"success": True, "info": "任务运行成功.", 'code':'C2300901'}
            task_info['content'] = {
                'ids': {
                    'id': "",
                    'name': ""
                    }}
            return json.dumps(task_info)
        else:
            task_info = {"success": False, "info": "任务运行失败.", 'code':'C2300902'}
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
        cmd += "i/metabolome/metabset_delete "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_type="1",
            # submit_location="ProteinsetVenn",
            proteinset_id="5a9de1715b8c915efb2a3293",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
