# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import glob
import unittest
import types
import json
from bson.objectid import ObjectId


class GeneFusionCircosWorkflow(Workflow):
    """
    Used gene fusion circos plot

    #chr_pos_file
    Chromosome      chromStart      chromEnd        Name    Stain
      2 chr1    0       2300000 p36.33  gneg
      3 chr1    2300000 5300000 p36.32  gpos25
      4 chr1    5300000 7100000 p36.31  gneg

    #gene_fusion_file
    left_gene       right_gene
    ENO1    GSTM1
    MPL     RPS15
    WT1     EGFR

    #gene_pos_file
    Chromosome      chromStart      chromEnd        Gene
      2 chr1    8921418 8934967 ENO1
      3 chr1    17345375        17380514        SDHB
      4 chr1    27022894        27107247        ARID1A
      5 chr1    41976121        42501596        HIVEP3
      6 chr1    43803519        43818443        MPL
      7 chr1    45794977        45805926        MUTYH
      8 chr1    65300244        65351947        JAK1
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GeneFusionCircosWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "chr_pos_file", "type": "infile", "format": "ref_rna_v2.common"},  # 记录染色体信息的文件
            {"name": "gene_pos_file", "type": "infile", "format": "ref_rna_v2.common"},  # 记录染色体信息的文件,
            {"name": "gene_fusion_file", "type": "infile", "format": "ref_rna_v2.common"},  # 记录染色体信息的文件,
            # {"name": "project_type", "type": "string","default": None},  # fasta文件
            {"name": "target_chrs", "type": "string", "default": ""},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.gene_fusion_circos")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(GeneFusionCircosWorkflow, self).run()


    def run_tool(self):
        opts = {
            'chr_pos_file': self.option("chr_pos_file"),
            'gene_pos_file': self.option("gene_pos_file"),
            'gene_fusion_file': self.option("gene_fusion_file"),
            'target_chrs': self.option("target_chrs"),
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        output_dir = self.tool.output_dir
        pdf_path = glob.glob(os.path.join(output_dir, '*.pdf'))
        pdf_s3_path = os.path.join(self._sheet.output,os.path.basename(pdf_path[0]))
        diff_ma = self.api.api("tool_lab.gene_fusion_circos")
        diff_ma.add_s3_result(
                                 main_id=self.option('main_id'),
                                 s3_output=pdf_s3_path,
                                 )
        self.set_output()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量差异ma图",0],
            [r'gene_fusion_circos.pdf', 'pdf', '基因融合圈图', 0],
            [r'diff_ma.png', 'png', '表达量差异ma图', 0],
            [r'diff_ma.svg', 'svg', '表达量差异ma图', 0],
        ])

        super(GeneFusionCircosWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.diff_ma import GeneFusionCircosWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "diff_volcano" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.diff_ma",
            "options": dict(
                raw_file='/mnt/ilustre/users/sanger-dev/workspace/20210322/DiffexpBatch_jgo0_k3guobgpi1lrvksgs4pu61_6989_3353/DiffexpBatch/output/ma.xls',
                pvalue=0.05,
                fc=2,
                x_axis_name="log10(TPM)",
                y_axis_name="log2(FC)",
                title_name="MA Plot",
                color="ref_blue_grey"
            )
        }
        wsheet = Sheet(data=data)
        wf =GeneFusionCircosWorkflow(wsheet)
        wf.sheet.id = 'diff_ma'
        wf.sheet.project_sn = 'diff_ma'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
