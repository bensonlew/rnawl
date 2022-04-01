# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class ExpNormWorkflow(Workflow):
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
        super(ExpNormWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "exp_matrix", "type": "infile", "format": "ref_rna_v2.common"},  # 表达量矩阵
            {"name": "gene_length", "type": "infile", "format": "ref_rna_v2.common"},  # 基因长度文件
            {'name': 'convert_type', 'type': 'string'},
            # {"name": "float_num", "type": "int", "default": '4'},  # 保留小数位数
            {'name': 'method_type', 'type': 'string'},
            {"name": "intersect", "type": "bool", "default": True},  # 是否取交集
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.exp_norm")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(ExpNormWorkflow, self).run()

    def check_options(self):
        self.convert = self.option('convert_type')
        if self.option('method_type') == 'in':
            if self.convert not in ["tpm", "cpm", "fpkm"]:
                raise OptionError('方法选择错误，请从tpm，cpm，fpkm中选取')
        if self.option('method_type') == 'out':
            if self.convert not in ["TMM", "TMMwzp", "RLF", "uqua", "DESeq2"]:
                raise OptionError('方法选择错误，请从TMM，RLF，uqua，DESeq2中选取')
        if self.convert not in ["tpm", "cpm", "fpkm", "TMM", "TMMwzp", "RLF", "uqua", "DESeq2"]:
            raise OptionError('不支持该标准化方法')
        if self.convert in ['tpm', 'fpkm']:
            if not self.option("gene_length").is_set:
                raise OptionError("{} should provie gene length file!".format(self.convert))
        return True

    def run_tool(self):

        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'convert_type': self.convert,
            # 'float_num': self.option('float_num'),
        }
        if self.convert in ['tpm', 'fpkm']:
            opts.update({'gene_length': self.option('gene_length'), 'intersect': self.option('intersect')})
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
            [".", "", "表达量标准化目录", 0],
        ])
        super(ExpNormWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.exp_norm import ExpNormWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "exp_norm" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.exp_norm",
            "options": dict(
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/known_seqs_count.matrix",
                convert_type="DESeq2",
                # float_num=4,
            )
        }
        wsheet = Sheet(data=data)
        wf =ExpNormWorkflow(wsheet)
        wf.sheet.id = 'exp_norm'
        wf.sheet.project_sn = 'exp_norm'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
