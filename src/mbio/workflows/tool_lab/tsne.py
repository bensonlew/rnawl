# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class TsneWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TsneWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'group', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.tsne")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(TsneWorkflow, self).run()

    def run_tool(self):
        opts = {
            'exp_file': self.option('exp_file').path,
            'group': self.option('group').path
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        # picedit = self.api.api('tool_lab.picedit')
        # picedit.add_translation(self.option('main_id'), self.tool.option('translation_file').path)
        s3_path = '{}/{}'.format(self._sheet.output, 'TSNE.pdf')
        tsne = self.api.api('tool_lab.tsne')
        tsne.add_tsne(self.option('main_id'), s3_path)
        self.set_output()


    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TSNE数据降维", 0],
        ])
        super(TsneWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.exp_norm import ExpNormWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "exp_norm" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.exp_norm",
            "options": dict(
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/known_seqs_count.matrix",
                convert_type="DESeq2",
                # float_num=4,
            )
        }
        wsheet = Sheet(data=data)
        wf =ExpNormWorkflow(wsheet)
        wf.sheet.id = 'exp_norm'
        wf.sheet.project_sn = 'exp_norm'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
