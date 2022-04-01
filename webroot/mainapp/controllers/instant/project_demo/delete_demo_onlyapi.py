# -*- coding: utf-8 -*-
import web
import json
import os
import unittest
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.core.base import Base
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo


class DeleteDemoAction(Base):
    def __init__(self):
        super(DeleteDemoAction, self).__init__()
        self.input_data = web.input()
        self._project_type = self.input_data.type

    @check_sig
    def POST(self):
        requires = ['type', 'task_id']
        for i in requires:
            if not (hasattr(self.input_data, i)):
                return json.dumps({"success": False, "info": "Lack argument: %s!" % i})

        find_result = self.db['sg_table_relation'].find_one({})
        if find_result:
            target = find_result['target']
            print(target)
        else:
            print("nothing found in sg_table_relation")
        try:
            del_demo = DeleteDemoMongo(self.input_data.task_id, self.input_data.type, target)
            del_demo.run()
            run_info = dict(success=True)
            run_info['info'] = 'delete_demo_of_{}'.format(self.input_data.task_id)
            return json.dumps(run_info)
        except Exception as e:
            return json.dumps({"success": False, "info": "running error: %s" % str(e)})


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/project_demo/delete_demo "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="copy_demo_ffff",
            type="ref_rna_v2"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
