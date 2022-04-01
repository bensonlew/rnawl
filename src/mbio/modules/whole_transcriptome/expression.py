# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.module import Module

from mbio.packages.whole_transcriptome.utils import read_fastq_dir


class ExpressionModule(Module):
    def __init__(self, work_id):
        super(ExpressionModule, self).__init__(work_id)
        PROGRAM = ('rsem', 'kallisto', 'salmon')
        options = [
            {'name': 'program', 'type': 'string', 'default': PROGRAM[2]},
            {'name': 'transcripts', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'strand_specific', 'type': 'bool', 'default': True},
            {'name': 'strand_dir', 'type': 'string', 'default': 'RF'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 't2g', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'filter', 'type': 'bool', 'default': True},
            {'name': 'threshold', 'type': 'float', 'default': 0.0},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(ExpressionModule, self).run()
        self.run_index()

    def run_index(self):
        self.index = self.add_tool('whole_transcriptome.expression.index')
        opts = {
            'program': self.option('program').lower(),
            'transcripts': self.option('transcripts'),
            't2g': self.option('t2g'),
        }
        self.index.set_options(opts)
        self.index.on('end', {'rsem': self.run_rsem, 'kallisto': self.run_kallisto, 'salmon': self.run_salmon}[
            self.option('program').lower()])
        self.index.run()

    def run_rsem(self):
        is_se, fastqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        reference_name = self.index.option('rsem_index').path
        for sample, fastqs in fastqs_dict.items():
            rsem = self.add_tool('whole_transcriptome.expression.rsem')
            opts = {'reference_name': reference_name, 'sample': sample}
            if self.option('strand_specific'):
                if self.option('strand_dir').startswith('R'):
                    opts['strandedness'] = 'reverse'
                elif self.option('strand_dir').startswith('F'):
                    opts['strandedness'] = 'forward'
            else:
                opts['strandedness'] = 'none'
            opts['upstream_read_file'] = fastqs[0]
            if not is_se:
                opts['downstream_read_file'] = fastqs[1]
            rsem.set_options(opts)
            self.tools.append(rsem)
        else:
            self.on_rely(self.tools, self.run_mashup)
        for tool in self.tools:
            tool.run()

    def run_kallisto(self):
        is_se, fastqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        index = self.index.option('kallisto_index').path
        for sample, fastqs in fastqs_dict.items():
            kallisto = self.add_tool('whole_transcriptome.expression.kallisto')
            opts = {'index': index, 'strand_specific': self.option('strand_specific'),
                    'strand_dir': self.option('strand_dir'), 't2g': self.option('t2g'), 'sample': sample}
            if is_se:
                opts.update({'single': True, 'reads': fastqs[0]})
            else:
                opts.update({'reads_1': fastqs[0], 'reads_2': fastqs[1]})
            kallisto.set_options(opts)
            self.tools.append(kallisto)
        else:
            self.on_rely(self.tools, self.run_mashup)
        for tool in self.tools:
            tool.run()

    def run_salmon(self):
        is_se, fastqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        index = self.index.option('salmon_index').path
        gene_map = self.option('t2g')
        for sample, fastqs in fastqs_dict.items():
            salmon = self.add_tool('whole_transcriptome.expression.salmon')
            opts = {'index': index, 'gene_map': gene_map, 'sample': sample}
            if is_se:
                opts.update({'unmated': fastqs[0]})
            else:
                if self.option('strand_specific'):
                    if self.option('strand_dir').startswith('R'):
                        lib_type = 'ISR'
                    elif self.option('strand_dir').startswith('F'):
                        lib_type = 'ISF'
                else:
                    lib_type = 'A'
                opts.update({'mates1': fastqs[0], 'mates2': fastqs[1], 'lib_type': lib_type})
            salmon.set_options(opts)
            self.tools.append(salmon)
        else:
            self.on_rely(self.tools, self.run_mashup)
        for tool in self.tools:
            tool.run()

    def run_mashup(self):
        self.mashup = self.add_tool('whole_transcriptome.expression.mashup')
        t_quant_list = os.path.join(self.work_dir, 'T.quant.list')
        g_quant_list = os.path.join(self.work_dir, 'G.quant.list')
        if self.option('program') == 'rsem':
            with open(t_quant_list, 'w') as t_handle, open(g_quant_list, 'w') as g_handle:
                for tool in self.tools:
                    t_handle.write('{}\t{}\n'.format(tool.option('sample'), tool.option('isoforms_results').path))
                    g_handle.write('{}\t{}\n'.format(tool.option('sample'), tool.option('genes_results').path))
        elif self.option('program') == 'kallisto':
            with open(t_quant_list, 'w') as t_handle, open(g_quant_list, 'w') as g_handle:
                for tool in self.tools:
                    t_handle.write('{}\t{}\n'.format(tool.option('sample'), tool.option('t_abundance').path))
                    g_handle.write('{}\t{}\n'.format(tool.option('sample'), tool.option('g_abundance').path))
        elif self.option('program') == 'salmon':
            with open(t_quant_list, 'w') as t_handle, open(g_quant_list, 'w') as g_handle:
                for tool in self.tools:
                    t_handle.write('{}\t{}\n'.format(tool.option('sample'), tool.option('quant_sf').path))
                    g_handle.write('{}\t{}\n'.format(tool.option('sample'), tool.option('quant_genes_sf').path))
        opts = {
            'program': self.option('program'),
            't_quant_list': t_quant_list,
            'g_quant_list': g_quant_list,
            'filter': self.option('filter'),
            'threshold': self.option('threshold'),
            't2g': self.option('t2g')
        }
        self.mashup.set_options(opts)
        self.mashup.on('end', self.set_output)
        self.mashup.run()

    def set_output(self):
        for fname in os.listdir(self.mashup.output_dir):
            source = os.path.join(self.mashup.output_dir, fname)
            link_name = os.path.join(self.output_dir, fname)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        self.end()

    def end(self):
        super(ExpressionModule, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.expression',
            'instant': False,
            'options': {
                'program': 'salmon',
                'transcripts': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/LargeGush/MergeKnownNew/transcript.fasta',
                'strand_specific': True,
                'strand_dir': 'RF',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/FastpRna/fastq_dir',
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/LargeGush/MergeKnownNew/t2g.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
