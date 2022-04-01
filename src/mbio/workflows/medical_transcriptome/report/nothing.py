# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,qinjincheng'

from biocluster.workflow import Workflow
import os
import shutil
import glob
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import unittest
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo

def workfuncdeco(func):
    def wrapper(*args, **kwargs):
        # print(str(len(*args)))
        args[0].logger.info('begin of the function ({}) at ({})'.format(func.__name__, func.__module__))
        result = func(*args, **kwargs)
        args[0].logger.info('final of the function ({}) at ({})'.format(func.__name__, func.__module__))
        # print(str(len(*args)))
        return result
    return wrapper

class NothingWorkflow(Workflow):
    '''
    last_modify: 2019.06.12
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NothingWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'opt1', 'type': 'string'},
            {'name': 'opt2', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self._sheet.output = self.option("s3_output")
        self.config.DBVersion = 1
        self.tool = self.add_tool("medical_transcriptome.snp.nothing")

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))


    def run(self):
        self.logger.info("DB_version ：{}".format(self.config.DBVersion))
        # self.tool.on("end", self.delete_mongo_data)
        self.start_listener()
        self.run_prepare()
        self.delete_mongo_data()

        # self.run_tool()
        # self.delete_mongo_data()
        super(NothingWorkflow, self).run()

    def delete_mongo_data(self):
        self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        self.program = os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin/python')
        a= DeleteDemoMongo("jiushiceshixia", 'ref_rna_v2')
        a.run()
        # cmd = '{} {}'.format(self.program, self.script)
        # cmd += ' {} {}'.format("jiushixiaceshi", 'ref_rna_v2')
        # code = os.system(cmd)
        # if code == 0:
        #     self.logger.info("命令{}执行成功！".format(cmd))
        # else:
        #     raise Exception("命令{}执行失败！".format(cmd))
        self.end()


    def run_prepare(self):
        self.logger.info("就给你们看看")


    def run_tool(self):
        opts={"opt1":self.option("opt1"),
              "opt2": self.option("opt2")}
        self.tool.set_options(opts)
        self.tool.run()

    def set_db(self):
        self.end()
        self.logger.info("看够了没有")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(NothingWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.nothing import NothingWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "nothing" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.report.nothing",
            "options": dict(
                opt1='oleilei',
                opt2='wulala',
            )
        }
        wsheet = Sheet(data=data)
        wf =NothingWorkflow(wsheet)
        wf.sheet.id = 'nothing'
        wf.sheet.project_sn = 'nothing'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

