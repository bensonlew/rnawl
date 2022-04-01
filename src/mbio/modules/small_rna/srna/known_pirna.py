# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
import sys
import os
import glob
import unittest
from Bio import SeqIO
import re


class KnownPirnaModule(Module):
    '''
    last_modify: 2019.04.11
    '''
    def __init__(self, work_id):
        super(KnownPirnaModule, self).__init__(work_id)
        options = [
            {'name': 'fa', 'type': 'infile', 'format': 'small_rna.fasta'},
            {'name': 'species', 'type': 'string'},
            {'name': 'config', 'type': 'infile', 'format': 'small_rna.common'}
        ]
        self.add_option(options)
        self.bowtie_tools = list()
        self.length_pirna = {'ssc': [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]}

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(KnownPirnaModule, self).run()
        self.run_tool()

    def run_tool(self):
        self.bowtie_pirna()
        for tool in self.bowtie_tools:
            tool.run()
        self.on_rely(self.bowtie_tools, self.stat)

    def fasta_split(self):
        fasta = self.option('fa').path
        pirna = SeqIO.parse(fasta, 'fasta')
        for record in pirna:
            length = len(str(record.seq))
            if length in self.length_pirna['ssc']:
                with open('{}/{}_len_{}.fa'.format(self.output_dir, self.option('species'), length), 'a+') as f:
                    f.write('>{}'.format(record.id) + '\n' + str(record.seq) + '\n')



    def bowtie_pirna(self):
        self.fasta_split()
        dir = self.output_dir
        print dir
        for home, dirs, files in os.walk(dir):
            for filename in files:
                fa_path = os.path.join(home, filename)
                length = re.findall(r"len_(.*).fa", fa_path)[0]
                self.bowtie_pirna = self.add_tool('small_rna.srna.known_pirna_new')
                options = {
                    'fa': fa_path,
                    'length': length,
                    'species': self.option('species'),
                    'config': self.option('config').path
                }
                self.bowtie_pirna.set_options(options)
                self.bowtie_tools.append(self.bowtie_pirna)


    def stat(self):
        self.stat = self.add_tool('small_rna.srna.pirna_stat')
        pirna_list = os.path.join(self.work_dir, 'pirna.list')
        with open(pirna_list, 'w') as handle:
            for module in self.bowtie_tools:
                handle.write('{}\n'.format(module.option('pirna_stat').path))
        options = {
            'pir_list': pirna_list,

        }
        self.stat.set_options(options)
        self.stat.on('end', self.set_output)
        self.stat.run()


    def set_output(self):
        self.end()

    def end(self):
        super(KnownPirnaModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bowtie_pirna_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'small_rna.srna.known_pirna',
            'instant': False,
            'options': {
                'fa': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test/uniq_10000.fasta',
                'species': 'ssc',
                'config': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test/qc_file.config'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    import sys
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
