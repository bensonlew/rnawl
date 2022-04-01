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


class GenesetDeleteAction(DenovoRnaV2Controller):
    def __init__(self):
        super(GenesetDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        #print data
        default_argu = ['geneset_id', "task_id"]
        for argu in default_argu:
            if not hasattr(data, argu):
                var = []
                var.append(argu)
                info = {'success': False, 'info': "Lack argument: %s" % (argu), "code": 'C1601101', "variables": var}
                return json.dumps(info)
        #print(data.geneset_id)
        if self.denovo_rna_v2.delete_geneset(data.geneset_id, data.task_id):
            task_info = {"success": True, "info": "任务运行成功."}
            task_info['content'] = {
                'ids': {
                    'id': "",
                    'name': ""
                    }}
            return json.dumps(task_info)
        else:
            task_info = {"success": False, "info": "任务运行失败.", "code": 'C1601102', "variables": ''}
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
        cmd += "i/denovo_rna_v2/geneset_delete "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="denovo_rna_v2",
            task_type="1",
            # submit_location="GenesetVenn",
            geneset_ids="5a35da8ba4e1af38e6fb83d7",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
