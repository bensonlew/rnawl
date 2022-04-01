# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import subprocess
import unittest

import pandas as pd
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class GetfastaAgent(Agent):
    '''
    last_modify: 2019.10.21
    '''

    def __init__(self, parent):
        super(GetfastaAgent, self).__init__(parent)
        options = [
            {'name': 'details', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},

            {'name': 'circbed', 'type': 'outfile', 'format': 'whole_transcriptome.common'},

            {'name': 'circfasta', 'type': 'outfile', 'format': 'whole_transcriptome.fasta'},

        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(GetfastaAgent, self).end()


class GetfastaTool(Tool):
    def __init__(self, config):
        super(GetfastaTool, self).__init__(config)
        self.program = {
            'bedtools': 'bioinfo/rna/bedtools2-master/bin/bedtools'

        }
        self.file = {
            'circbed': os.path.join(self.output_dir, 'circrna.bed'),
            'circfasta': os.path.join(self.output_dir, 'circrna.fasta'),

        }
        self.script = {

        }

    def run(self):
        super(GetfastaTool, self).run()
        self.run_bed()
        self.run_getfasta()
        self.fasta()
        self.set_output()
        self.end()

    def run_bed(self):
        df_bed = dict()
        df = pd.read_table(self.option('details').path)
        df_bed['chr'] = df['chr']
        df_bed['start'] = df['circrna_start']
        df_bed['end'] = df['circrna_end']
        df_bed['name'] = df['circrna_id']
        df_bed['score'] = 0
        df_bed['strand'] = df['strand']
        circ_bed = pd.DataFrame(df_bed)
        circ_bed_reindex = circ_bed.reindex(columns=['chr', 'start', 'end', 'name', 'score', 'strand'])
        circ_bed_reindex.to_csv(self.file['circbed'], index=False, header=False, sep='\t')

    def run_getfasta(self):
        if os.path.isfile(self.file['circfasta']):
            n_records = len(list(SeqIO.parse(self.file['circfasta'], 'fasta')))
            n_bedlines = len(open(self.file['circbed']).readlines())
            if n_records == n_bedlines:
                return
        cmd = '{} getfasta -fi {} -bed {} -fo {} -name -s'.format(self.program['bedtools'], self.option('genome').path,
                                                         self.file['circbed'], self.file['circfasta'])
        cmd_name = 'run_getfasta'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code == -11:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    def fasta(self):
        cmd = "sed -i 's/(.*)//g' {}".format(self.file['circfasta'])
        self.logger.info(cmd)
        os.system(cmd)
        self.logger.info("done")

        # cmd_name = 'run_fasta'
        # command = self.add_command(cmd_name, cmd, shell=True)
        # command.software_dir = 'user/bin'
        # command.run()
        # self.wait()
        # if command.return_code == 0:
        #     self.logger.info("{} 运行成功".format(cmd_name))
        # elif command.return_code is None:
        #     self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
        #     command.rerun()
        #     self.wait()
        #     if command.return_code is 0:
        #         self.logger.info("{} 运行成功".format(cmd_name))
        #     else:
        #         self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        # else:
        #     self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

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
            'name': 'whole_transcriptome.circrna.getfasta',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Gallus_gallus/GRCg6a_ensembl_101_/dna/Gallus_gallus.dna.toplevel.fa',
                'details': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/longRNA-seq/detail.txt'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
