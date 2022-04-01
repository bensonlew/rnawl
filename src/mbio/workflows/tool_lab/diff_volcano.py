# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import glob


class DiffVolcanoWorkflow(Workflow):
    """
    Used for cds to protein code
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffVolcanoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "raw_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            {"name": "pvalue", "type": "float", "default": 0.05},
            {"name": "fc", "type": "float", "default": 2.0},
            #翻译方向 ["5","3","both"]
            {"name": "x_axis_name", "type": "string", "default": "log2(FC)"},
            # 翻译起始位点 ["1","2","3"]
            {"name": "y_axis_name", "type": "string", "default": "-log10(Padjust)"},
            {"name": "title_name", "type": "string","default": "Volcano Plot"},
            {"name": "color", "type": "string", "default": "red_blue_grey"},
            {"name": "showname", "type": "string", "default": "no"},
            {"name": "showname_num", "type": "string", "default": ""},
            {"name": "showname_str", "type": "string", "default": ""},
            {"name": "main_id", "type": "string",}, # 序列转化的方法，默认c，即直接取["c","r","cr"]
            # {"name": "min_len", "type": "int", },
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.diff_volcano")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(DiffVolcanoWorkflow, self).run()

    def run_tool(self):
        opts = {
            'raw_file': self.option("raw_file").prop["path"],
            'pvalue': self.option('pvalue'),
            'fc': self.option('fc'),
            'x_axis_name' : self.option("x_axis_name"),
            'y_axis_name': self.option("y_axis_name"),
            'title_name': self.option("title_name"),
            'color': self.option("color"),
        }
        if self.option("showname") == "yes":
            opts.update({
                "showname_num" : self.option("showname_num"),
                "showname_str" : self.option("showname_str"),
                })
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        output_dir = self.tool.output_dir
        pdf_path = glob.glob(os.path.join(output_dir, '*.pdf'))
        pdf_s3_path = os.path.join(self._sheet.output,os.path.basename(pdf_path[0]))
        diff_volcano = self.api.api("tool_lab.diff_volcano")
        diff_volcano.add_s3_result(
                                 main_id=self.option('main_id'),
                                 s3_output=pdf_s3_path,
                                 )
        self.set_output()
        # pdf_list = list()
        # for i in pdf_path:
        #     pdf_list.append('{}/{}'.format(self._sheet.output, os.path.basename(i)))


    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量差异火山图",0],
            [r'.*\.volcano\.pdf', 'xls', '表达量差异火山图', 0],
        ])
        super(DiffVolcanoWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.diff_volcano import DiffVolcanoWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "diff_volcano" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.diff_volcano",
            "options": dict(
                raw_file='/mnt/ilustre/users/sanger-dev/workspace/20210322/DiffexpBatch_jgo0_k3guobgpi1lrvksgs4pu61_6989_3353/DiffexpBatch/output/voca.xls',
                pvalue=0.05,
                fc=2,
                x_axis_name="log2(FC)",
                y_axis_name="-log10(Padjust)",
                title_name="Volcano Plot",
                color="ref_blue_grey"
            )
        }
        wsheet = Sheet(data=data)
        wf =DiffVolcanoWorkflow(wsheet)
        wf.sheet.id = 'diff_volcano'
        wf.sheet.project_sn = 'diff_volcano'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
