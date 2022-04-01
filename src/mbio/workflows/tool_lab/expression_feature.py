# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class ExpressionFeatureWorkflow(Workflow):

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpressionFeatureWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'sam', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'bamorsam', 'type': 'string', 'default': 'sam'},
            {'name': 'gtforgfforsaf', 'type': 'string', 'default': 'gtf'},
            {'name': 'gtf', 'type': 'string', 'format': 'whole_transcriptome.common'},
            {'name': 'allowMultiOverlap', 'type': 'bool', 'default': False},
            {'name': 'countMultiMappingReads', 'type': 'bool', 'default': False},
            {'name': 'fraction', 'type': 'bool', 'default': False},
            {'name': 'isPairedEnd', 'type': 'bool', 'default': True},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.feature_count")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(ExpressionFeatureWorkflow, self).run()

    def check_options(self):
        # self.convert = self.option('convert_type')
        # if self.convert not in ["tpm", "cpm", "fpkm", "TMM", "TMMwzp", "RLF", "uqua", "DESeq2"]:
        #     raise OptionError('不支持该标准化方法')
        # if self.convert in ['tpm', 'fpkm']:
        #     if not self.option("gene_length").is_set:
        #         raise OptionError("{} should provie gene length file!".format(self.convert))
        # return True
        pass
    def run_tool(self):
        opts = {
            'sam': self.option('sam'),
            'gtf': self.option('gtf'),
            'allowMultiOverlap': self.option('allowMultiOverlap'),
            'countMultiMappingReads': self.option('countMultiMappingReads'),
            'fraction': self.option('fraction'),
            'isPairedEnd': self.option('isPairedEnd')
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量", 0],
        ])
        super(ExpressionFeatureWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.expression_feature import ExpressionFeatureWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "expression_feature" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.expression_feature",
            "options":   {
                'sam': '/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/HTSeq-0.6.1/scripts/file.sam',
                'gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf',
                }
        }
        wsheet = Sheet(data=data)
        wf =ExpressionFeatureWorkflow(wsheet)
        wf.sheet.id = 'expression_feature'
        wf.sheet.project_sn = 'expression_feature'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
