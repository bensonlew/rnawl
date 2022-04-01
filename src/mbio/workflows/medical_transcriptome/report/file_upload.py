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


class FileUploadWorkflow(Workflow):
    '''
    last_modify: 2019.06.12
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FileUploadWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 's3_output', 'type': 'string'},
            {'name': 'file_path', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self.option("s3_output")

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def run(self):
        self.start_listener()
        self.run_tool()
        super(FileUploadWorkflow, self).run()

    def run_tool(self):
        if os.path.isdir(self.option("file_path")):
            for file in os.listdir(self.option("file_path")):
                if os.path.isfile(os.path.join(self.option("file_path"),file)):
                    raw_file = os.path.join(self.option("file_path"),file)
                    new_file = os.path.join(self.output_dir, file)
                    os.link(raw_file ,new_file )
                elif os.path.isdir(os.path.join(self.option("file_path"),file)):
                    raw_file = os.path.join(self.option("file_path"), file)
                    new_file = os.path.join(self.output_dir, file)
                    shutil.copytree(raw_file,new_file)
        elif os.path.isfile(self.option("file_path")):
            raw_file = self.option("file_path")
            new_file = os.path.join(self.output_dir,os.path.basename(self.option("file_path")))
            os.link(raw_file, new_file)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(FileUploadWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.file_upload import FileUploadWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "file_upload" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.report.file_upload",
            "options": dict(
                # s3_output='s3://medicaltranscriptome/files/m_188/45psjlpdvgn6qulm2he6oepva2/6sieeiseu0ubhjnsm1o6302l5r/intermediate_results/SequenceDetail/',
                # file_path='/mnt/ilustre/users/sanger-dev/workspace/20201207/Single_rmats_4546_9201/Detail/output/detail/pep_seq',
                s3_output='s3://medicaltranscriptome/files/m_188/45psjlpdvgn6qulm2he6oepva2/emvhefm6875rraot9k7c0jkv70/workflow_results/03Gene_structure_analysis/01AS/',
                file_path='/mnt/ilustre/users/sanger-dev/workspace/20201222/MedicalTranscriptome_emvhefm6875rraot9k7c0jkv70/upload/03Gene_structure_analysis/01AS'
            )
        }
        wsheet = Sheet(data=data)
        wf =FileUploadWorkflow(wsheet)
        wf.sheet.id = 'file_upload'
        wf.sheet.project_sn = 'file_upload'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

