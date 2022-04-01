# -*- coding: utf-8 -*-
# __author__ = 'shijin'
import web
import json
import os
import unittest
import datetime
from mainapp.models.mongo.denovo_rna_v2 import DenovoRnaV2
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController


class ProteinsetDeleteAction(ItraqTmtController):
    def __init__(self):
        super(ProteinsetDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        #print data
        default_argu = ['proteinset_id', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': "Lack argument: %s", "variables":[argu], "code" : "C1901001"}
                return json.dumps(info)
        #print(data.geneset_id)
        if self.itraq_tmt.delete_proteinset(data.proteinset_id, data.task_id):
            task_info = {"success": True, "info": "success"}
            task_info['content'] = {
                'ids': {
                    'id': "",
                    'name': ""
                    }}
            return json.dumps(task_info)
        else:
            task_info = {"success": False, "info": "failed", "code" : "C1901002"}
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
        cmd += "i/itraq_and_tmt/proteinset_delete "
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
