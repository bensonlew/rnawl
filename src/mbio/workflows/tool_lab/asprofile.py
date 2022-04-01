# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class AsprofileWorkflow(Workflow):
    """
    可变剪切事件

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AsprofileWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gtf_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_options(self):
        # for k, v in self.sheet.options().items():
        #     self.logger.debug('{} = {}'.format(k, v))
        pass

    def run(self):
        self.run_asprofile()
        super(AsprofileWorkflow, self).run()

    def run_asprofile(self):
        self.asprofile = self.add_module('tool_lab.asprofile')
        self.asprofile.set_options({
            'gtf_dir': self.option('gtf_dir'),
            'ref_fa': self.option('ref_fa')
        })
        self.asprofile.on('end', self.set_db)
        self.asprofile.run()

    def end(self):
        super(AsprofileWorkflow, self).end()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        asprofile = self.api.api("tool_lab.asprofile")
        # add result info
        as_result = os.path.join(self.asprofile.output_dir, 'AS_result_merge.txt')
        as_statistics = os.path.join(self.asprofile.output_dir, 'AS_statistics_merge.txt')
        main_id = asprofile.add_asprofile_result(as_result, self.option('main_id'))
        asprofile.add_asprofile_statistics(as_statistics, main_id)
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.asprofile import AsprofileWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "ASprofile" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.asprofile",
            "options": dict(
                gtf_dir='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/gtf/',
                ref_fa='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/ref.fa'
            )
        }

        wsheet = Sheet(data=data)
        wf =AsprofileWorkflow(wsheet)
        wf.sheet.id = 'ASprofile'
        wf.sheet.project_sn = 'ASprofile'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
