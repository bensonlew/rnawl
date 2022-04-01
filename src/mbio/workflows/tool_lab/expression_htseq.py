# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class ExpressionHtseqWorkflow(Workflow):

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpressionHtseqWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'bamorsam', 'type': 'string', 'default': 'bam'},
            {'name': 'sam', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'gtforgff', 'type': 'string', 'default': 'gtf'},
            {'name': 'gtf', 'type': 'string', 'format': 'whole_transcriptome.common'},
            {'name': 'method', 'type': 'string', 'default': 'union'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.htseq_count")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(ExpressionHtseqWorkflow, self).run()

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
            'bamorsam': self.option('bamorsam'),
            'gtf': self.option('gtf'),
            'method': self.option('method')
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
        super(ExpressionHtseqWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.expression_htseq import ExpressionHtseqWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "expression_htseq" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.expression_htseq",
            "options":   {
                'sam': '/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/HTSeq-0.6.1/scripts/file.sam',
                'gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf',
                'method': "union",
                'bamorsam': 'sam'
                }
        }
        wsheet = Sheet(data=data)
        wf =ExpressionHtseqWorkflow(wsheet)
        wf.sheet.id = 'expression_rsem'
        wf.sheet.project_sn = 'expression_rsem'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
