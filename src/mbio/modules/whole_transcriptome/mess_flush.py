# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.module import Module


class MessFlushModule(Module):
    def __init__(self, work_id):
        super(MessFlushModule, self).__init__(work_id)
        PROGRAM = ('rsem', 'kallisto', 'salmon')
        options = [
            {'name': 'program', 'type': 'string', 'default': PROGRAM[2]},
            {'name': 'transcripts', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'strand_specific', 'type': 'bool', 'default': True},
            {'name': 'strand_dir', 'type': 'string', 'default': 'RF'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 't2g', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 't_style', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_style', 'type': 'infile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(MessFlushModule, self).run()
        self.run_expression()

    def run_expression(self):
        self.expression = self.add_module('whole_transcriptome.expression')
        transcripts = self.option('transcripts').path
        t2g = self.option('t2g').path
        opts = {
            'program': self.option('program'),
            'transcripts': transcripts,
            'strand_specific': self.option('strand_specific'),
            'strand_dir': self.option('strand_dir'),
            'fastq_dir': self.option('fastq_dir'),
            't2g': t2g
        }
        self.expression.set_options(opts)
        self.expression.on('end', self.run_exp_style)
        self.expression.run()

    def run_exp_style(self):
        self.exp_style = self.add_tool('whole_transcriptome.formation.exp_style')
        opts = {
            't_tpm': os.path.join(self.expression.output_dir, 'T.tpm.txt'),
            'g_tpm': os.path.join(self.expression.output_dir, 'G.tpm.txt'),
            't_fpkm': os.path.join(self.expression.output_dir, 'T.fpkm.txt'),
            'g_fpkm': os.path.join(self.expression.output_dir, 'G.fpkm.txt'),
            't_style': self.option('t_style'),
            'g_style': self.option('g_style')
        }
        self.exp_style.set_options(opts)
        self.exp_style.on('end', self.set_output, 'exp_style')
        self.exp_style.run()

    def set_output(self, event):
        obj = event['bind_object']
        if os.path.isdir(os.path.join(self.output_dir, event['data'])):
            shutil.rmtree(os.path.join(self.output_dir, event['data']))
        shutil.copytree(obj.output_dir, os.path.join(self.output_dir, event['data']))
        self.end()

    def end(self):
        super(MessFlushModule, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mess_flush_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.mess_flush',
            'instant': False,
            'options': {
                'program': 'kallisto',
                'transcripts': '/mnt/lustre/users/sanger/workspace/20200106/Single_assembly_9312_9082/Assembly/output/all.fasta',
                'strand_specific': True,
                'strand_dir': 'RF',
                'fastq_dir': '/mnt/lustre/users/sanger/workspace/20191230/Longrna_sanger_231556/FastpRna/output/fastq',
                't2g': '/mnt/lustre/users/sanger/workspace/20200106/Single_assembly_9312_9082/Assembly/output/t2g.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
