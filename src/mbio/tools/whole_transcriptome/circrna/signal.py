# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd
from Bio import SeqIO

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class SignalAgent(Agent):
    '''
    last_modify: 2019.8.22
    '''

    def __init__(self, parent):
        super(SignalAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'circRNA_merge', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'bed', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'fasta', 'type': 'outfile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'signal', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(SignalAgent, self).end()


class SignalTool(Tool):
    def __init__(self, config):
        super(SignalTool, self).__init__(config)
        self.program = {
            'bedtools': 'bioinfo/rna/bedtools2-master/bin/bedtools',
        }
        self.file = {
            'fasta': os.path.join(self.output_dir, '{}_circRNA.fasta'.format(self.option('sample'))),
            'bed': os.path.join(self.output_dir, '{}_circRNA.bed'.format(self.option('sample'))),
            'signal': os.path.join(self.output_dir, '{}_circ_merge_signal.txt'.format(self.option('sample')))
        }


    def run(self):
        super(SignalTool, self).run()
        self.run_make_bed()
        self.run_bedtools()
        self.run_signal()
        self.set_output()
        self.end()

    def run_make_bed(self):
        f = pd.read_table(self.option('circRNA_merge').path)
        a = dict()
        a['chr'] = f['chr']
        a['circRNA_start'] = f['circRNA_start'] - 3
        a['circRNA_end'] = f['circRNA_end'] + 2
        bed = pd.DataFrame(a)
        bed = bed[['chr', 'circRNA_start', 'circRNA_end']]
        bed.to_csv(self.file['bed'],index = False,sep = '\t',header = False)

    def run_bedtools(self):
        cmd = '{} getfasta -fi {} -bed {} -fo {}'.format(self.program['bedtools'], self.option('genome').path,self.file['bed'], self.file['fasta'])
        cmd_name = 'run_bedtools'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            pass
        elif command.return_code == -11:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running run_bedtools")
                pass
        else:
            self.set_error("fail to run run_bedtools")

        # runcmd(self,cmd_name,cmd)

    def run_signal(self):
        signal = []
        id = []
        for seq_records in SeqIO.parse(self.file['fasta'], 'fasta'):
            a = str(seq_records.seq[-2:]) + str(seq_records.seq[0:2])
            signal.append(a)
            id.append(seq_records.id)
        file = pd.read_table(self.option('circRNA_merge').path)
        file['Junction_site_type'] = signal
        file.to_csv(self.file['signal'], index = False, sep = '\t')

    def set_output(self):
        self.option('signal').set_path(self.file['signal'])



class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'signal_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.signal',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'circRNA_merge': '/mnt/ilustre/users/sanger-dev/workspace/20190924/Single_merge_6217_4019/Merge/output/circRNA_merge.txt'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
