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


class GeneFusionAnnotWorkflow(Workflow):
    """
    Used for cds to protein code
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GeneFusionAnnotWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "raw_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            {"name": "project_type", "type": "string", "default": None},  # fasta文件
            # {"name": "raw_file_info", "type": "string"},  # fasta文件
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.gene_fusion_annot")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(GeneFusionAnnotWorkflow, self).run()

    def get_raw_file_info(self):
        if not self.option("project_type"):
            raw_file_path = self.option("raw_file").prop["path"]
            project_type = "custom"
            return project_type, raw_file_path
        elif  self.option("project_type") == "medical" :
            raw_file_path = self.option("raw_file").prop["path"]
            project_type = "medical_transcriptome"
            return project_type, raw_file_path

    # def get_raw_file_info(self):
    #     try:
    #         raw_file_info = json.loads(self.option("raw_file_info"))
    #         project_type = raw_file_info["project_type"]
    #         raw_file_path = raw_file_info["file_path"]
    #         return project_type,raw_file_path
    #     except:
    #         project_type = "custom"
    #         raw_file_path = self.option("raw_file_info")
    #         return project_type,raw_file_path



    def run_tool(self):
        project_type,file_path = self.get_raw_file_info()
        opts = {
            'raw_file': file_path,
            'project_type' :project_type,
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        output_dir = self.tool.output_dir
        result_file_path = os.path.join(output_dir,"Fusion_anno.xls")
        fusion_annot = self.api.api("tool_lab.gene_fusion_annot")
        fusion_annot.add_fusion_annot(
                                 main_id=self.option('main_id'),
                                 result_file_path = result_file_path,
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
            [r'Fusion_anno.xls', 'xls', '基因融合注释分析结果表', 0],
        ])
        super(GeneFusionAnnotWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.diff_ma import GeneFusionAnnotWorkflow
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
        wf =GeneFusionAnnotWorkflow(wsheet)
        wf.sheet.id = 'diff_ma'
        wf.sheet.project_sn = 'diff_ma'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
