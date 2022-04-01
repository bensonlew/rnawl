# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class ExtractGffFastaWorkflow(Workflow):
    """
    用于提取gff/gtf文件中的mRNA、cds以及pep序列
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExtractGffFastaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "genome", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因文件
            {"name": "gff", "type": "infile", "format": "ref_rna_v2.gtf"},  # 参考基因的注释文件
            {"name": "strand", "type": "string", "default": ''},  # + or -
            {"name": "chr", "type": "string", "default": ''},  # 指定染色体名称
            {"name": "start", "type": "int", "default": ''},  # 起始位置
            {"name": "end", "type": "int", "default": ''},  # 终止为止
            {"name": "custom_position", "type": "bool", "default": False},  # 提取特定染色体指定位置序列
            {"name": "seq_type", "type": "string", "default": ''},  # 序列类型, mRNA/CDS/protein/all
            {"name": "reverse", "type": "bool", "default": False},
            {"name": "no_pseudo", "type": "bool", "default": False},  # 过滤掉含有 'pseudo' 的注释信息
            {"name": "no_cds", "type": "bool", "default": False},  # 丢弃掉无CDS的转录本
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.extract_gff_fasta")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(ExtractGffFastaWorkflow, self).run()

    def run_tool(self):
        opts = {
            'genome': self.option('genome'),
            'gff': self.option('gff'),
            'seq_type': self.option('seq_type'),
            'no_pseudo': self.option('no_pseudo'),
            'no_cds': self.option('no_cds'),
            'reverse': self.option('reverse'),
        }
        if self.option("custom_position"):
            opts.update({'chr': self.option('chr'),
                         'strand': self.option('strand'),
                         'start': self.option('start'),
                         'end': self.option('end'),
                         })
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
            [".", "", "提取基因组序列文件",0],
            ["mRNA.fa", "", "转录本序列",0],
            ["cds.fa", "", "cds序列", 0],
            ["protein.fa", "", "蛋白序列",0],
        ])
        super(ExtractGffFastaWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.extract_gff_fasta import ExtractGffFastaWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "extract_gff_fasta" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.extract_gff_fasta",
            "options": dict(
                genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                gff=test_dir + "/" + "gtf/Saccharomyces_cerevisiae.R64-1-1.39.gtf",
                seq_type='all',
                reverse=False,
                no_cds=False,
                chr="avs",
                strand="+",
                start=10,
                end=10000
            )
        }
        wsheet = Sheet(data=data)
        wf =ExtractGffFastaWorkflow(wsheet)
        wf.sheet.id = 'extract_gff_fasta'
        wf.sheet.project_sn = 'extract_gff_fasta'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
