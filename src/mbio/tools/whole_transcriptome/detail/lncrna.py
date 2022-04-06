# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.tool import Tool
import subprocess
from biocluster.core.exceptions import FileError
import time
from mbio.packages.whole_transcriptome.utils import runcmd
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout



class LncrnaAgent(Agent):
    '''
    last_modify: 2019.10.31
    '''

    def __init__(self, parent):
        super(LncrnaAgent, self).__init__(parent)
        options = [
            {'name': 'lncrna_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'known_lncrna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'novel_lncrna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'dna_fa', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'transcript', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'relation_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'biomart_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'biomart_type', 'type': 'string', 'default': None},
            {'name': 'ref_cds', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'ref_pep', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new_cds', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new_pep', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'genome_id', 'type': 'string', 'default': None},
            {'name': 'transcript_bed', 'type': 'outfile', 'format': 'gene_structure.bed'},
            {'name': 'transcript_fa', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'gene_bed', 'type': 'outfile', 'format': 'gene_structure.bed'},
            {'name': 'gene_fa', 'type': 'outfile', 'format': 'sequence.fasta'}
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('dna_fa').path) / 1024 ** 3 + 16)

    def end(self):
        super(LncrnaAgent, self).end()


class LncrnaTool(Tool):
    def __init__(self, config):
        super(LncrnaTool, self).__init__(config)
        self.program = {
            'gffread': 'bioinfo/rna/cufflinks-2.2.1/gffread',
            'gtftogenepred': 'bioinfo/align/ucsc_tools/gtfToGenePred',
            'genepredtobed': 'bioinfo/align/ucsc_tools/genePredToBed',
            'python': 'miniconda2/bin/python',
            'bedtools': 'bioinfo/ref_rna_v2/miniconda2/bin/bedtools'
        }
        self.script = {
            'get_gene_bed': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/get_gene_bed.py'),
            'fasta_clean': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/detail/fasta_clean.py'),
            'base_build': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/detail/lnc_base_build.py')
        }
        self.file = {
            'transcript_fasta': os.path.join(self.output_dir, 'transcript.fasta'),
            'transcript_genepred': os.path.join(self.work_dir, 'transcript.genepred'),
            'transcript_bed': os.path.join(self.output_dir, 'transcript.bed'),
            'gene_bed': os.path.join(self.output_dir, 'gene.bed'),
            'gene_fasta': os.path.join(self.output_dir, 'gene.fasta'),
            'input_json': os.path.join(self.work_dir, 'input.json'),
            'database': os.path.join(self.output_dir, 'lncrna.db'),
            'gene_detail': os.path.join(self.output_dir, 'gene_detail.pk'),
            'transcript_detail': os.path.join(self.output_dir, 'transcript_detail.pk')
        }
        self.gtf2bed_path = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/gtf2bed.pl')
        self.gtf2standard = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/gtf_standard.pl')
        self.perl = 'program/perl-5.24.0/bin/perl'
        self.process_rerun = 0

    def run(self):
        super(LncrnaTool, self).run()
        self.run_gffread()
        self.run_gtftogenepred()
        self.run_genepredtobed()
        # self.run_gtftobed()
        self.run_get_gene_bed()
        if self.option("transcript").is_set:
            self.file['transcript_fasta'] = self.option("transcript").path
            pass
        else:
            self.run_bedtools_getfasta()
        self.run_fasta_clean()
        self.pre_base_build()
        self.run_base_build()
        self.set_output()
        self.end()

    def run_gffread(self):
        cmd = '{} {} -T'.format(self.program['gffread'], self.option('lncrna_gtf').path)
        cmd += ' -g {}'.format(self.option('dna_fa').path)
        cmd += ' -w {}'.format(self.file['transcript_fasta'])
        runcmd(self, 'run_gffread', cmd, block=False)

    def run_gtftogenepred(self):
        cmd = '{} {} {} -ignoreGroupsWithoutExons'.format(
            self.program['gtftogenepred'], self.option('lncrna_gtf').path, self.file['transcript_genepred']
        )
        runcmd(self, 'run_gtftogenepred', cmd)

    def run_genepredtobed(self):
        cmd = '{} {} {}'.format(self.program['genepredtobed'], self.file['transcript_genepred'],
                                self.file['transcript_bed'])
        runcmd(self, 'run_genepredtobed', cmd)

    def run_gtftobed(self):
        cmd = "perl {} --i {} >  {}".format(self.gtf2standard , self.option('lncrna_gtf').path, self.option('lncrna_gtf').path + ".std")
        try:
            subprocess.check_output(cmd, shell=True)
            print cmd
        except subprocess.CalledProcessError:
            raise FileError("Operation error")
        pass
        # runcmd(self, 'run_gtfstandard', cmd)

        cmd = "{} {} {} {}".format(self.perl, self.gtf2bed_path, self.option('lncrna_gtf').path + ".std", self.file['transcript_bed'])
        runcmd(self, 'run_gtftobed', cmd)

    def run_get_gene_bed(self):
        cmd = '{} {}'.format(self.program['python'], self.script['get_gene_bed'])
        cmd += ' --gtf {}'.format(self.option('lncrna_gtf').path)
        cmd += ' --bed {}'.format(self.file['transcript_bed'])
        cmd += ' --output {}'.format(self.file['gene_bed'])
        runcmd(self, 'run_get_gene_bed', cmd)

    def run_bedtools_getfasta(self):
        cmd = '{} getfasta -name -s'.format(self.program['bedtools'])
        cmd += ' -fi {}'.format(self.option('dna_fa').path)
        cmd += ' -fo {}'.format(self.file['gene_fasta'])
        cmd += ' -bed {}'.format(self.file['gene_bed'])
        runcmd(self, 'run_bedtools_getfasta', cmd)

    def run_fasta_clean(self):
        cmd = '{} {}'.format(self.program['python'], self.script['fasta_clean'])
        cmd += ' -i {}'.format(self.file['gene_fasta'])
        cmd += ' -o {}'.format(self.file['gene_fasta'])
        runcmd(self, 'run_fasta_clean', cmd)

    def pre_base_build(self):
        database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
        try:
            collection = database['sg_genome_db']
            genome_doc = collection.find_one({'genome_id': self.option('genome_id')})
            organism_name = genome_doc['organism_name']
            species_urls = genome_doc['ensemble_web']
            if 'ensembl' in genome_doc['ensemble_web'].lower():
                source = 'ensembl'
            elif 'ncbi' in genome_doc['ensemble_web'].lower():
                source = 'ncbi'
            else:
                source = None
            json.dump({'lncrna_gtf': self.option('lncrna_gtf').path,
                       'biomart_file': self.option('biomart_file').path,
                       'biomart_type': self.option('biomart_type'),
                       'relation_file': self.option('relation_file').path,
                       'known_lncrna_detail': self.option('known_lncrna_detail').path,
                       'novel_lncrna_detail': self.option('novel_lncrna_detail').path,
                       'transcript_fa': self.file['transcript_fasta'],
                       'transcript_bed': self.file['transcript_bed'],
                       'organism_name': organism_name,
                       'species_urls': species_urls,
                       'source': source,
                       'transcript_detail_pk': self.file['transcript_detail']},
                      open(self.file['input_json'], 'w'), indent=4)
        except (ServerSelectionTimeoutError, NetworkTimeout): # 捕获因为mongo服务器问题导致的异常后重运行此方法
            if self.process_rerun < 5:
                self.process_rerun += 1
                self.logger.info("检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun))
                time.sleep(5)
                self.pre_base_build()
            else:
                self.add_state('memory_limit', '检测到TimeoutError, 重运行tool')
    def run_base_build(self):
        cmd = '{} {}'.format(self.program['python'], self.script['base_build'])
        cmd += ' --json {}'.format(self.file['input_json'])
        cmd += ' --database {}'.format(self.file['database'])
        runcmd(self, 'run_base_build', cmd)

    def set_output(self):
        self.option('transcript_bed').set_path(self.file['gene_fasta'])
        self.option('transcript_fa').set_path(self.file['transcript_fasta'])
        self.option('gene_bed').set_path(self.file['gene_bed'])
        self.option('gene_fa').set_path(self.file['gene_fasta'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'lncrna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.detail.lncrna',
            'instant': False,
            'options': {
                'lncrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/filter_by_express/filtered_file/all_lncrna.gtf',
                'dna_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'biomart_file': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/biomart/biomart.txt',
                'biomart_type': 'type1',
                'genome_id': 'GM0259',
                'known_lncrna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/known_lnc_identify/known_lncrna_detail.xls',
                'novel_lncrna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/new_lncrna_predict/novel_lncrna_predict_detail.xls',
                'relation_file': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/filter_by_express/filtered_file/trans_type.xls'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
