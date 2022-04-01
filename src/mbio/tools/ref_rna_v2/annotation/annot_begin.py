# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest

class AnnotBeginAgent(Agent):
    '''
    last_modify: 2019.04.22
    '''
    def __init__(self, parent):
        super(AnnotBeginAgent, self).__init__(parent)
        options = [
            {'name': 'gene2trans', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'gtf', 'type': 'infile', 'format': 'ref_rna_v2.gtf'},
            {'name': 'fasta', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'g2t2p', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 't2g2r2l2p', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'longest_t2g', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps('annot_begin')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.annot_begin.start()
        self.step.update()

    def step_end(self):
        self.step.annot_begin.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(AnnotBeginAgent, self).end()

class AnnotBeginTool(Tool):
    def __init__(self, config):
        super(AnnotBeginTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.gffread = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/cufflinks-2.2.1/gffread')
        self.get_representative_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/get_representative.py')
        self.get_relation_map_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/get_relation_map.py')
        self.t2g2r2l2p_tsv = os.path.join(self.output_dir, 't2g2r2l2p.tsv')
        self.longest_t2g_tsv = os.path.join(self.output_dir, 'longest.t2g.tsv')

    def run(self):
        super(AnnotBeginTool, self).run()
        if self.option('gene2trans').is_set:
            self.run_get_representative()
        elif self.option('gtf').is_set and self.option('fasta').is_set and self.option('g2t2p').is_set:
            self.run_get_relation_map()
        self.set_output()
        self.end()

    def run_get_representative(self):
        shutil.copy(self.option('gene2trans').prop['path'], self.t2g2r2l2p_tsv)
        cmd = '{} {}'.format(self.python, self.get_representative_py)
        cmd += ' --input {}'.format(self.t2g2r2l2p_tsv)
        cmd += ' --output {}'.format(self.longest_t2g_tsv)
        cmd_name = 'run_get_representative'
        self.run_code(cmd_name, cmd)

    def run_get_relation_map(self):
        cmd = '{} {}'.format(self.python, self.get_relation_map_py)
        cmd += ' --gffread {}'.format(self.gffread)
        cmd += ' --gtf {}'.format(self.option('gtf').path)
        cmd += ' --fasta {}'.format(self.option('fasta').path)
        cmd += ' --g2t2p {}'.format(self.option('g2t2p').path)
        cmd += ' --output {}'.format(self.output_dir)
        cmd_name = 'run_get_relation_map'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run %s, abord', variables=(n), code="33710102")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('t2g2r2l2p').set_path(self.t2g2r2l2p_tsv)
        self.option('longest_t2g').set_path(self.longest_t2g_tsv)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_begin_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.annotation.annot_begin',
            'instant': False,
            'options': {
                'gtf': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'fasta': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
                'g2t2p': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/Annotation_v2/g2t2p',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
