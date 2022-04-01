# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
import os
import unittest
import pandas as pd
from Bio import SeqIO
import re

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class GetfastaCircAgent(Agent):
    '''
    last_modify: 2019.10.21
    '''

    def __init__(self, parent):
        super(GetfastaCircAgent, self).__init__(parent)
        options = [
            {'name': 'circbase', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'circbase_bed', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},



            {'name': 'circfasta', 'type': 'outfile', 'format': 'ref_rna_v2.fasta'},


        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '40G'

    def end(self):
        super(GetfastaCircAgent, self).end()


class GetfastaCircTool(Tool):
    def __init__(self, config):
        super(GetfastaCircTool, self).__init__(config)
        self.program = {
            'bedtools': 'bioinfo/rna/bedtools2-master/bin/bedtools'


        }
        self.file = {
            'circbase_bed':os.path.join(self.output_dir,'circbase_bed.bed'),
            'circfasta': os.path.join(self.output_dir, 'circrna.fasta'),

        }
        self.script = {


        }


    def run(self):
        super(GetfastaCircTool, self).run()
        self.run_circbase_bed()
        # self.run_bed()
        self.run_getfasta()
        self.set_output()
        self.end()

    def run_circbase_bed(self):
        with open(self.file['circbase_bed'], 'w') as bed:
            for seq_record in SeqIO.parse(self.option('circbase').path, "fasta"):
                id = seq_record.id
                (circbase_id, chr_start_end_strant, host_gene_id, protein_id) = id.rstrip().split('|')
                chr_start_end = chr_start_end_strant[:-1]
                (chr, start, end) = re.split('[:-]', chr_start_end)
                bed.write(chr+'\t'+str(int(start)+1)+'\t'+end+'\n')




    #
    #
    #
    # def run_bed(self):
    #     df_bed = dict()
    #     df = pd.read_table(self.option('circbase').path)
    #     df_bed['chr']=df['chr']
    #     df_bed['start']=df['circrna_start']
    #     df_bed['end']=df['circrna_end']
    #     circ_bed = pd.DataFrame(df_bed)
    #     circ_bed_reindex = circ_bed.reindex(columns=['chr', 'start', 'end'])
    #     circ_bed_reindex.to_csv(self.file['circbed'], index=False, header=False, sep='\t')

    def run_getfasta(self):
        cmd = '{} getfasta -fi {} -bed {} -fo {}'.format(self.program['bedtools'], self.option('genome').path, self.file['circbase_bed'], self.file['circfasta'])
        cmd_name = 'run_getfasta'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('circfasta').set_path(self.file['circfasta'])





class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'getfasta_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.getfasta_circ',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/database/mouse/Mus_mm9.fasta',
                'circbase': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/mouse_mm9_circRNAs_putative_spliced_sequence.fa'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


