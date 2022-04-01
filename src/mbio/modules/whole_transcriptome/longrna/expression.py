# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import unittest

from biocluster.module import Module


class ExpressionModule(Module):
    def __init__(self, work_id):
        super(ExpressionModule, self).__init__(work_id)
        options = [
            {'name': 'program', 'type': 'string', 'default': 'salmon'},
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
            'program': self.option('program'),
            'transcripts': self.option('transcripts'),
        }
        self.index.set_options(opts)
        self.index.on('end', self.run_salmon)
        self.index.run()

    def run_salmon(self):
        is_se, sample2fastq_list_dict = read_fastq_dir(self.option('fastq_dir').path)
        index = self.index.option('salmon_index').path
        gene_map = self.option('t2g')
        for sample, fastq_list in sample2fastq_list_dict.items():
            salmon = self.add_tool('whole_transcriptome.expression.salmon')
            opts = {'index': index, 'gene_map': gene_map, 'sample': sample}
            if is_se:
                opts.update({'unmated': fastq_list[0]})
            else:
                if self.option('strand_specific'):
                    if self.option('strand_dir').startswith('R'):
                        lib_type = 'ISR'
                    elif self.option('strand_dir').startswith('F'):
                        lib_type = 'ISF'
                else:
                    lib_type = 'A'
                opts.update({'mates1': fastq_list[0], 'mates2': fastq_list[1], 'lib_type': lib_type})
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
        for file_name in os.listdir(self.mashup.output_dir):
            source = os.path.join(self.mashup.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        else:
            self.end()

    def end(self):
        super(ExpressionModule, self).end()


def read_fastq_dir(fastq_dir):
    is_se = False
    sample2fastq_list_dict = dict()
    for line in open(os.path.join(fastq_dir, 'list.txt')):
        str_list = line.strip().split('\t')
        fastq = os.path.join(fastq_dir, str_list[0])
        sample, mate_type = str_list[1:]
        if sample in sample2fastq_list_dict:
            if mate_type == 'l':
                sample2fastq_list_dict[sample].insert(0, fastq)
            elif mate_type == 'r':
                sample2fastq_list_dict[sample].append(fastq)
        else:
            sample2fastq_list_dict[sample] = [fastq]
    else:
        pe_sample_count = len(filter(lambda item: len(item[1]) > 1, sample2fastq_list_dict.items()))
        if pe_sample_count:
            if pe_sample_count == len(sample2fastq_list_dict):
                is_se = False
            else:
                raise Exception('mix mate type found in {} -> {}'.format(fastq_dir, sample2fastq_list_dict))
        else:
            is_se = True
        return is_se, sample2fastq_list_dict


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
                'transcripts': '/mnt/ilustre/users/sanger-dev/workspace/20190906/Single_large_gush_5347_9539/LargeGush/MergeKnownNew/output/known_and_new.fa',
                'strand_specific': True,
                'strand_dir': 'RF',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/clean_data/little',
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20190906/Single_large_gush_5347_9539/LargeGush/MergeKnownNew/output/t2g.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
