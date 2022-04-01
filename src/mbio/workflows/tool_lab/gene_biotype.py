# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
import re
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from biocluster.config import Config


class GeneBiotypeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GeneBiotypeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input_type", "type": "string", "default": ""},  # gtf or gff
            {"name": "gtf", "type": "infile", "format": "small_rna.common"},  # 输入gtf/gff文件
            {"name": "gene_list", "type": "infile", "format": "small_rna.common"},  # 目标基因集列表
            {"name": "biotype", "type": "string", "default": ""},  # biotype类型
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.gene_biotype")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(GeneBiotypeWorkflow, self).run()

    def check_options(self):
        for type in self.option('biotype').split(","):
            if type.lower() not in ['all', 'mrna', 'mirna', 'lncrna', 'trna', 'rrna', 'snrna', 'pseudogene']:
                raise OptionError('{} not supported right row.'.format(type))
        if self.option('input_type') not in ['gff', 'gtf']:
            raise OptionError('input file type {} is not supported.'.format(self.option('input_type')))
        return True

    def run_tool(self):
        options = {
            "input_type": self.option("input_type"),
            "gtf": self.option("gtf"),
            "gene_list": self.option("gene_list"),
            "biotype": self.option("biotype"),
        }
        self.tool.set_options(options)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Gene_Biotype提取结果目录", 0],
            ["gene_biotype.xls", "xls", "基因Biotype提取结果文件", 0],
        ])
        super(GeneBiotypeWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.gene_biotype import GeneBiotypeWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            "id": "GeneBiotype_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.gene_biotype",
            "options": dict(
                input_type='gff',
                gtf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/tmp/gene_biotype/GCF_000001405.39_GRCh38.p13_genomic_partial.gff",
                biotype="lncRNA,miRNA",
                gene_list="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/tmp/gene_biotype/gene_list"
                # input_type='gtf',
                # gtf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/tmp/gene_biotype/GCF_000001405.39_GRCh38.p13_genomic_partial.gtf",
                # biotype="All",
            )
        }
        wsheet = Sheet(data=data)
        wf =GeneBiotypeWorkflow(wsheet)
        wf.sheet.id = 'gene_biotype'
        wf.sheet.project_sn = 'gene_biotype'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
