# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.config import Config


class StrWorkflow(Workflow):
    """
    短串重复序列

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(StrWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'fq_type', 'type': "string", "default": "PE"},
            {'name': 'length_required', 'type': 'string', 'default': '30'},
            {'name': 'quality_score_system', 'type': 'string', 'default': 'phred 33'},
            {'name': 'adapter_a', 'type': 'string', 'default': 'AGATCGGAAGAGCACACGTC'},
            {'name': 'adapter_b', 'type': 'string', 'default': 'AGATCGGAAGAGCGTCGTGT'},
            {'name': 'genome_id', 'type': 'string'},
            {'name': 'strand_specific', 'type': 'bool', 'default': False},
            {'name': 'ref_genome', 'type': 'string', 'default': None},
            {'name': 'sample_list', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'variant-catalog', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
        collection = database['sg_genome_db']
        genome_info = collection.find_one({'genome_id': self.option('genome_id')})
        self.genome_version = genome_info['assembly']
        self.annot_version = genome_info['annot_version']
        self.dna_fa = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish/{}'.format(genome_info["dna_fa"]))
        self.tools = list()
        self.samples = list()


    def check_options(self):
        # for k, v in self.sheet.options().items():
        #     self.logger.debug('{} = {}'.format(k, v))
        pass

    def run(self):
        self.run_fastp_rna()
        super(StrWorkflow, self).run()

    def run_fastp_rna(self):
        self.fastp_rna = self.add_module('tool_lab.fastp_rna')
        fq_dir = self.option('fastq_dir').path
        sample_path = os.path.join(fq_dir, 'abs.list.txt')
        open(sample_path, 'w').writelines(
            '{}/{}'.format(fq_dir, line) for line in open(os.path.join(fq_dir, 'list.txt'))
        )
        if self.option('fq_type') == "PE":
            self.fastp_rna.set_options({
                'sample_path': sample_path,
                'fq_type': self.option('fq_type'),
                'length_required': self.option('length_required'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_sequence': self.option('adapter_a'),
                'adapter_sequence_r2': self.option('adapter_b')
            })
        else:
            self.fastp_rna.set_options({
                'sample_path': sample_path,
                'fq_type': self.option('fq_type'),
                'length_required': self.option('length_required'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_sequence_s': self.option('adapter_a'),
            })
        self.fastp_rna.on('end', self.run_rnaseq_mapping)
        self.fastp_rna.run()

    def run_rnaseq_mapping(self):
        fastq_dir = self.fastp_rna.option('sickle_dir')
        self.rnaseq_mapping = self.add_module('ref_rna_v2.rnaseq_mapping')
        self.rnaseq_mapping.set_options({
            'ref_genome': self.option('ref_genome'),
            'genome_version': self.genome_version,
            'genome_annot_version': self.annot_version,
            'mapping_method': 'hisat',
            'seq_method': self.option('fq_type'),
            'fastq_dir': fastq_dir,
            'assemble_method': 'stringtie',
            'strand_specific': self.option('strand_specific')
        })
        self.rnaseq_mapping.on('end', self.run_str_expansion_hunter)
        self.rnaseq_mapping.run()

    def run_str_expansion_hunter(self):
        bamlist = os.path.join(self.rnaseq_mapping.output_dir, 'bamlist')
        self.str_expansion = self.add_module('tool_lab.str_expansion')
        self.str_expansion.set_options({
            'bamlist': bamlist,
            'ref': self.option('ref_genome'),
            'variant-catalog': self.option('variant-catalog')
        })
        self.str_expansion.on('end', self.end)
        self.str_expansion.run()



    def run_str(self):
        bamlist = os.path.join(self.rnaseq_mapping.output_dir, 'bamlist')
        self.str_predict = self.add_module('tool_lab.str')
        self.str_predict.set_options({
            'bamlist': bamlist,
            'ref_fa': self.dna_fa,
            'sample_list': self.option('sample_list')
        })
        self.str_predict.on('end', self.set_db)

    def end(self):
        super(StrWorkflow, self).end()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # asprofile = self.api.api("tool_lab.asprofile")
        # # add result info
        # as_result = os.path.join(self.asprofile.output_dir, 'AS_result_merge.txt')
        # as_statistics = os.path.join(self.asprofile.output_dir, 'AS_statistics_merge.txt')
        # main_id = asprofile.add_asprofile_result(as_result, self.option('main_id'))
        # asprofile.add_asprofile_statistics(as_statistics, main_id)
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.str import StrWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "STR" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.asprofile",
            "options": dict(
                fastq_dir='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/data/test_data_1',
                genome_id='GM0502',
                ref_genome='Solanum_lycopersicum',
                sample_list='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/data/sample_list.txt'
            )
        }

        wsheet = Sheet(data=data)
        wf =StrWorkflow(wsheet)
        wf.sheet.id = 'Str'
        wf.sheet.project_sn = 'Str'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
