# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
import zipfile
import glob

class ExpressionSalmonWorkflow(Workflow):

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpressionSalmonWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'fa_input', 'type': 'infile', 'format':'ref_rna_v2.fasta'},
            {'name': 'fq_type', 'type': 'string', 'default': 'PE'},
            {'name': 'fq_dir', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'fq_list', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 't2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'kmer', 'type': 'int', 'default': 31},
            {'name': 'phred_quals', 'type': 'int', 'default': 33},
            {'name': 'expression_table', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.module = self.add_module("tool_lab.expression")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(ExpressionSalmonWorkflow, self).run()

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
        fq_dir = os.path.dirname(self.option('fq_dir').path)
        print fq_dir
        fq_zip = glob.glob(fq_dir + '/*zip')[0]
        fz = zipfile.ZipFile(fq_zip, 'r')
        for f in fz.namelist():
            fz.extract(f, fq_dir)
        fq_new_list = os.path.join(self.work_dir, 'fq_new.list')
        with open(self.option('fq_list').path, 'r') as fq_list, open(fq_new_list, 'w') as fq_new:
            for line in fq_list.readlines():
                if self.option('fq_type') == 'SE':
                    sample, fq1 = line.strip().split('\t')
                    fq_new.write(sample + '\t' + os.path.join(fq_dir, fq1) + '\n')
                if self.option('fq_type') == 'PE':
                    sample, fq1, fq2 = line.strip().split('\t')
                    fq_new.write(sample + '\t' + os.path.join(fq_dir, fq1) + '\t' + os.path.join(fq_dir, fq2) + '\n')
        opts = {
            'fq_type': self.option('fq_type'),
            'fq_list': fq_new_list,
            'fa_input': self.option('fa_input'),
            'txpt2gene': self.option('t2g'),
            'method': 'salmon',
            'kmer': self.option('kmer')
        }
        self.module.set_options(opts)
        self.module.on('end', self.set_output)
        self.module.run()

    def set_output(self):
        for file in os.listdir(self.module.output_dir):
            os.link(os.path.join(self.module.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量", 0],
        ])
        super(ExpressionSalmonWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.expression_salmon import ExpressionSalmonWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "expression_salmon" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.expression_salmon",
            "options":   {'fq_type': 'PE',
                'fq_list': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/expression/fq.list',
                'fq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/expression/clean_data',
                'fa_input': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/transcript.fa',
                't2g': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-stringtie/output/NewTranscripts/trans2gene',
                'phred_quals': 33
                }
        }
        wsheet = Sheet(data=data)
        wf =ExpressionSalmonWorkflow(wsheet)
        wf.sheet.id = 'expression_salmon'
        wf.sheet.project_sn = 'expression_salmon'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
