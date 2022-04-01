# -*- coding: utf-8 -*-
# __author__ = 'shijin'
import web
import json
import datetime
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os

class GenesetDeleteAction(ProkRNAController):
    def __init__(self):
        super(GenesetDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        #print data
        default_argu = ['geneset_id', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        #print(data.geneset_id)
        if self.prok_rna.delete_geneset(data.geneset_id, data.task_id):
            task_info = {"success": True, "info": "任务运行成功."}
            task_info['content'] = {
                'ids': {
                    'id': "",
                    'name': ""
                    }}
            return json.dumps(task_info)
        else:
            task_info = {"success": False, "info": "任务运行失败."}
            return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/prok_rna/geneset_delete "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
               task_id="tsg_35470",
               geneset_id="5e14483d17b2bf244a1e7c97",
             )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()