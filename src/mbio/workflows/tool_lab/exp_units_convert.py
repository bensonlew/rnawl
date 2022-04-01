# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class ExpUnitsConvertWorkflow(Workflow):
    """
    Used for convert RNA-seq expression units convert, count/CPM/RPM/TPM/FPKM/RPKM are implemented.

    RPM/CPM: Reads/Counts of exon model per million mapped reads
    RPM/CPM=Total exon reads/ Mapped reads(Millions)

    RPKM/FPKM: Reads/Fragments Per Kilobase of exon model per Million mapped reads
    RPKM/FPKM=Total exon reads/[Mapped reads(Millions)*Exon length(Kb)]

    TPM is like RPKM and FPKM, except the order of operation is switched.
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpUnitsConvertWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "exp_matrix", "type": "infile", "format": "ref_rna_v2.common"},  # 表达量矩阵
            {"name": "gene_length", "type": "infile", "format": "ref_rna_v2.common"},  # 基因长度文件
            {"name": "from_unit", "type": "string", "default": ''},
            {"name": "to_unit", "type": "string", "default": ''},
            {"name": "float_num", "type": "int", "default": '4'},  # 保留小数位数
            {"name": "intersect", "type": "bool", "default": True},  # 是否取交集
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.exp_units_convert")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(ExpUnitsConvertWorkflow, self).run()

    def check_options(self):
        self.convert = self.option("from_unit").lower() + "2" + self.option("to_unit").lower()
        if self.convert not in ['count2cpm', 'fpkm2tpm', 'count2tpm', 'count2fpkm', 'cpm2fpkm', 'cpm2tpm', 'fpkm2cpm']:
            raise OptionError("%s to %s not supported!" % (self.option("from_unit").lower(), self.option("to_unit").lower()))
        if self.convert not in ['count2cpm', 'fpkm2tpm']:
            if not self.option("gene_length").is_set:
                raise OptionError("%s to %s should provie gene length file!") % (self.option("from_unit").lower(), self.option("to_unit").lower())
        return True

    def run_tool(self):

        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'gene_length': self.option('gene_length'),
            'convert_type': self.convert,
            'float_num': self.option('float_num'),
            'intersect': self.option('intersect'),
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
            [".", "", "基因表达量定量指标转换结果文件",0],
            [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
            [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
            [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
            [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
            [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
            [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
            [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        ])
        super(ExpUnitsConvertWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.extract_gff_fasta import ExpUnitsConvertWorkflow
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
