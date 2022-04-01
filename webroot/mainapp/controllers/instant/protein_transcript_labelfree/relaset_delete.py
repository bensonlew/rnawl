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
from mainapp.controllers.project.labelfree_controller import LabelfreeController


class RelasetDeleteAction(LabelfreeController):
    def __init__(self):
        super(RelasetDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        #print data
        default_argu = ['relaset_id', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        #print(data.geneset_id)
        if self.labelfree.delete_relaset(data.relaset_id, data.task_id):
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
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/protein_transcript_labelfree/relaset_delete "
        cmd += "-b http://bcl.tsg.com:9090 "
        args = dict(
            task_id='tsg_33097',
            task_type="1",
            # submit_location="ProteinsetVenn",
            relaset_id="5c734aacd887f30c07000029",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
