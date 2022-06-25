# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd
from mbio.packages.ref_rna.trans_step import merged_add_code


class TransBuildAgent(Agent):
    '''
    last_modify: 2019.10.08
    '''

    def __init__(self, parent):
        super(TransBuildAgent, self).__init__(parent)
        options = [
            {'name': 'program', 'type': 'string', 'default': 'gffcompare'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},
            {'name': 'seq_path', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'verbose', 'type': 'bool', 'default': True},
            {'name': 'merged_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(os.path.getsize(self.option('seq_path').path) / 1024.0 ** 3 * 2 + 20))

    def end(self):
        super(TransBuildAgent, self).end()


class TransBuildTool(Tool):
    def __init__(self, config):
        super(TransBuildTool, self).__init__(config)
        self.program = {
            'perl': 'miniconda2/bin/perl',
            'gffcompare': 'bioinfo/rna/gffcompare-0.9.8.Linux_x86_64/gffcompare',
            'cuffcompare': 'bioinfo/rna/cufflinks-2.2.1/cuffcompare',
            'python': 'miniconda2/bin/python',
            'gffread': 'bioinfo/rna/cufflinks-2.2.1/gffread',
        }
        self.script = {
            'gtfmerge': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/scripts/gtfmerge.pl'),
            'assembly_stat': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/scripts/assembly_stat.py'),
            'feature_class': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/assembly/feature_class.py'),
            'step_code': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/assembly/step_code.py')
        }
        self.file = {
            'add_code_merged_gtf': os.path.join(self.work_dir, 'add_code_merged.gtf'),
            'change_id_merged_gtf': os.path.join(self.work_dir, 'change_id_merged.gtf'),
            'ref_and_new_gtf': os.path.join(self.work_dir, 'ref_and_new.gtf'),
            'merged_gtf': os.path.join(self.work_dir, 'merged.gtf'),
            'tmap': os.path.join(self.work_dir, 'gffcmp.merged.gtf.tmap'),
            'ref_gtf': os.path.join(self.output_dir, 'ref.gtf'),
            'new_gtf': os.path.join(self.output_dir, 'new.gtf'),
            'all_gtf': os.path.join(self.output_dir, 'all.gtf'),
            'ref_fasta': os.path.join(self.output_dir, 'ref.fasta'),
            'new_fasta': os.path.join(self.output_dir, 'new.fasta'),
            'all_fasta': os.path.join(self.output_dir, 'all.fasta'),
        }

    def run(self):
        super(TransBuildTool, self).run()
        if self.option('program') == 'gffcompare':
            self.run_gffcompare()
            self.cal_merged_add_code()
            self.change_id()
            self.get_changed_new_gtf()
            self.concat_ref_and_new()
        self.run_feature_class()
        self.run_gffread_ref()
        self.run_gffread_new()
        self.run_gffread_all()
        self.run_step_code()
        self.set_output()
        self.end()

    def cal_merged_add_code(self):
        merged_add_code(trans_file=self.option('merged_gtf').path, tmap_file=self.file['tmap'],
                        new_trans=self.file['add_code_merged_gtf'])

    def change_id(self):
        cmd = '{} {}'.format(self.program['perl'], self.script['gtfmerge'])
        cmd += ' -i {}'.format(self.file['add_code_merged_gtf'])
        cmd += ' -compare {}'.format(self.file['tmap'])
        cmd += ' -ref {}'.format(self.option('ref_gtf').path)
        cmd += ' -o {}'.format(self.file['change_id_merged_gtf'])
        runcmd(self, 'change_id', cmd)

    def get_changed_new_gtf(self):
        if not os.path.isdir(os.path.join(self.work_dir, 'changed')):
            os.mkdir(os.path.join(self.work_dir, 'changed'))
        cmd = '{} {}'.format(self.program['python'], self.script['assembly_stat'])
        cmd += ' -tmapfile {}'.format(self.file['tmap'])
        cmd += ' -transcript_file {}'.format(self.file['change_id_merged_gtf'])
        cmd += ' -out_new_trans {}'.format(os.path.join(self.work_dir, 'changed', 'new_trans.gtf'))
        cmd += ' -out_new_genes {}'.format(os.path.join(self.work_dir, 'changed', 'new_genes.gtf'))
        cmd += ' -out_old_trans {}'.format(os.path.join(self.work_dir, 'changed', 'old_trans.gtf'))
        cmd += ' -out_old_genes {}'.format(os.path.join(self.work_dir, 'changed', 'old_genes.gtf'))
        runcmd(self, 'get_changed_new_gtf', cmd)

    def concat_ref_and_new(self):
        with open(self.file['ref_and_new_gtf'], 'w') as handle:
            handle.write(open(self.option('ref_gtf').path).read())
            handle.write(open(os.path.join(self.work_dir, 'changed', 'new_trans.gtf')).read())

    def run_gffcompare(self):
        if os.path.isfile(self.file['merged_gtf']):
            os.remove(self.file['merged_gtf'])
        os.link(self.option('merged_gtf').path, self.file['merged_gtf'])
        cmd = '{} '.format(self.program['gffcompare'])
        cmd += ' -r {}'.format(self.option('ref_gtf').path)
        cmd += ' -s {}'.format(self.option('seq_path').path)
        if self.option('verbose'):
            cmd += ' -V'
        cmd += ' {}'.format(self.file['merged_gtf'])
        runcmd(self, 'run_gffcompare', cmd)

    def run_feature_class(self):
        cmd = '{} {}'.format(self.program['python'], self.script['feature_class'])
        cmd += ' -r {}'.format(self.option('ref_gtf').path)
        cmd += ' -m {}'.format(self.file['ref_and_new_gtf'])
        cmd += ' -t {}'.format(self.file['tmap'])
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_feature_class', cmd)

    def run_gffread_ref(self):
        cmd = '{} {} -T'.format(self.program['gffread'], self.file['ref_gtf'])
        cmd += ' -g {}'.format(self.option('seq_path').path)
        cmd += ' -w {}'.format(self.file['ref_fasta'])
        runcmd(self, 'run_gffread_ref', cmd)

    def run_gffread_new(self):
        cmd = '{} {} -T'.format(self.program['gffread'], self.file['new_gtf'])
        cmd += ' -g {}'.format(self.option('seq_path').path)
        cmd += ' -w {}'.format(self.file['new_fasta'])
        runcmd(self, 'run_gffread_new', cmd)

    def run_gffread_all(self):
        cmd = '{} {} -T'.format(self.program['gffread'], self.file['all_gtf'])
        cmd += ' -g {}'.format(self.option('seq_path').path)
        cmd += ' -w {}'.format(self.file['all_fasta'])
        runcmd(self, 'run_gffread_all', cmd)

    def run_step_code(self):
        cmd = '{} {}'.format(self.program['python'], self.script['step_code'])
        cmd += ' -f {}'.format(self.file['all_fasta'])
        cmd += ' -t {}'.format(self.file['tmap'])
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_step_code', cmd)

    def run_cuffcompare(self):
        cmd = '{} -C'.format(self.program['cuffcompare'])
        cmd += ' -r {}'.format(self.option('ref_gtf').path)
        cmd += ' -s {}'.format(self.option('seq_path').path)
        if self.option('verbose'):
            cmd += ' -V'
        cmd += ' {}'.format(self.option('merged_gtf').path)
        runcmd(self, 'run_cuffcompare', cmd)

    def set_output(self):
        pass


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'trans_build_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.assembly.trans_build',
            'instant': False,
            'options': {
                'ref_gtf': '/mnt/ilustre/users/isanger/workspace/20200819/Longrna_majorbio_278257/FileCheck/haemonchus_contortus.PRJEB506.WBPS11.gtf',
                'seq_path': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/metazoa/Haemonchus_contortus/WormBase/dna/haemonchus_contortus.PRJEB506.WBPS11.genomic.fa',
                'merged_gtf': '/mnt/ilustre/users/isanger/workspace/20200819/Longrna_majorbio_278257/Assembly/TransBuild/merged.gtf',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
