# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import gevent
from bson.objectid import ObjectId
from mbio.packages.tool_lab.file_compress.file_compress import FileCompress

class SvgConvertWorkflow(Workflow):
    """
    Used for conver svg to jpg/png/tif/gif

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SvgConvertWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "svg_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            {"name": "target_type", "type": "string", "default": "png"},  # 目标格式，默认png ["png","tif","gif","jpg"
            {"name": "dpi", "type": "int", "default": "300"},  # 分辨率，默认300
            # {"name": "set_weight", "type": "bool", "default": False},  # 是否设置图片大小
            # {"name": "weight", "type": "int", "default": "0"},  # 改变图片宽度
            {"name": "sample_num", "type": "string", "default": "single"},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.svg_convert")
        self.set_options(self._sheet.options())
        self.tools=[]
        # self.weight=self.option("weight") if self.option("set_weight") else 0

    def run(self):
        self.run_tool()
        super(SvgConvertWorkflow, self).run()

    def run_tool(self):
        if self.option("sample_num") == "single":
            opts = {
                'svg_file': self.option('svg_file').prop["path"],
                'target_type': self.option('target_type'),
                'dpi': self.option('dpi'),
                # 'weight': self.weight,
                # 'min_len': self.option('min_len')
            }
            self.tool.set_options(opts)
            self.tool.on('end', self.set_output)
            self.tool.run()
        else:
            os.mkdir(os.path.join(self.work_dir,"svg_files"))
            unpress=FileCompress()
            file_name = os.path.basename(self.option('svg_file').prop["path"])
            unpress.uncompress(self.option('svg_file').prop["path"],os.path.join(self.work_dir,"svg_files"))
            # os.system("tar -xzvf {} -C {}".format(self.option("svg_file").prop["path"],os.path.join(self.work_dir,"svg_files")))
            svg_files=os.listdir(os.path.join(self.work_dir,"svg_files"))
            for file in sorted(svg_files):
                svg_convert=self.add_tool("tool_lab.svg_convert")
                svg_convert.set_options({
                    "svg_file": os.path.join(os.path.join(self.work_dir,"svg_files"),file),
                    'target_type': self.option('target_type'),
                    'dpi': self.option('dpi'),
                    # 'weight': self.weight,
                })
                self.tools.append(svg_convert)
            if self.tools:
                if len(self.tools) > 1:
                    self.on_rely(self.tools, self.set_output)
                elif len(self.tools) == 1:
                    self.tools[0].on('end', self.set_output)
                for tool in self.tools:
                    gevent.sleep(1)
                    tool.run()

    def set_output(self):
        if self.option("sample_num") == "single":
            for file in os.listdir(self.tool.output_dir):
                os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        else:
            for tool in self.tools:
                for file in os.listdir(tool.output_dir):
                    os.link(os.path.join(tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "svg转图片结果文件",0],
            # [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
            # [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
            # [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
            # [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
            # [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
            # [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
            # [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        ])
        super(SvgConvertWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.svg_convert import SvgConvertWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "svg_convert" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.svg_convert",
            "options": dict(
                # genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                # fasta_seq=">scaffold1\nATTATAAAATGAGGAGATGTAAATTTTAAGGAAAAAAATAAAGCTGTCGAGATTTTCTCGACAGCCTGAAGAGCACCCATCATAAAAGGGTGCTCTTTTTCTTTTCTCTTCCCGATTTGTTCACAGGCTGAATCCTCTCCTCATATGCTGAAAGGGAGATTCAGAAAATATTACAGACTTTTTTAAATAAAGAAGGTGAACGGTCAGTATGTTGAGCGGTTTAACGGTTGCGGTGATCGGGGGAGATGCAAGGCAGCTTGAAATCATTCGCAAGCTGTCACAGCAGCATGCCAAAGTGTTTTTGGTCGGATTTGATCAGCTGGATCATGGGTTTATCGGTGCTGAAAAGCTTAAAATGTCAGAACTTCCATTTGAACAGGTAGACAGTATGATTCTGCCGGTATCAGGTGCAACAGATGAAGGCGTCGTCGCCACAGTTTTCTCAAATGAGCAGGTCGTGCTGGAAGCAGAATATTTAGAAAGAACTCCAGCACATTGTACCTTGTACTCAGGTATTTCTAATACGTACTTAGACAATCTGGCAAAGCAGGTGAACCGGAAGCTTGTGAAGCTGTTTGAGCGCGATGATATTGCCATATATAACTCTATTCCAACAGTTGAAGGGATTATCATGATGGCCATTCAGCAAACGGACTATACGATTCATGGATCACATGTCGCTGTCCTCGGGCTTGGGAGAACAGGGCTCACAATTGCCCGCACAT",
                # svg_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/svg/input/chart.svg",
                # fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_filter/input/example.fasta",
                # fasta_file="",
                # set_weight=True,
                sample_num="multiple",
                svg_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/svg/input/test/test.tar.gz",
                target_type="png",
                dpi=350,
                # weight=200
                # max_len='10000'
            )
        }
        wsheet = Sheet(data=data)
        wf =SvgConvertWorkflow(wsheet)
        wf.sheet.id = 'extract_gff_fasta'
        wf.sheet.project_sn = 'extract_gff_fasta'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
