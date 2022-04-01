# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import re
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class RpmAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(RpmAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'rnasam', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name':'BSJ', 'type':'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'flagstat', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name':'RPM', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(RpmAgent, self).end()


class RpmTool(Tool):
    def __init__(self, config):
        super(RpmTool, self).__init__(config)
        self.program = {
            'samtools': 'bioinfo/align/samtools-1.3.1/samtools',
        }


        self.file = {
            'RPM': os.path.join(self.output_dir, '{}_circRNA_merge_signal_type_bsj_rpm.txt'.format(self.option('sample'))),
            'flagstat': os.path.join(self.output_dir,'{}_flagstat.txt'.format(self.option('sample')))

        }

    def run(self):
        super(RpmTool, self).run()
        self.run_samtools()
        self.run_rpm()
        self.set_output()
        self.end()

    def run_samtools(self):
        cmd = '{} flagstat {} > {}'.format(self.program['samtools'],self.option('rnasam').path,self.file['flagstat'])
        cmd_name = 'run_samtools'
        runcmd(self,cmd_name,cmd,shell=True)

    def run_rpm(self):
        with open(self.file['flagstat']) as f:
            for lines in f.readlines():
                result = re.findall(r"(.*)\s\+\s\d\smapped",lines)
                if result:
                    mappingreads = int(result[0])
                    print mappingreads

        circ = pd.read_table(self.option('BSJ').path)
        RPM = list()
        for i in range(0,len(circ)):
            rpm = float(1000000*circ['junction_reads'][i])/float((sum(circ['junction_reads']))+mappingreads)
            RPM.append(rpm)

        circ['RPM'] = RPM
        circ.to_csv(self.file['RPM'], index=False,sep = '\t')


    def set_output(self):
        self.option('RPM').set_path(self.file['RPM'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rpm_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.rpm',
            'instant': False,
            'options': {
                'sample': 'zjx',
                'rnasam': '/mnt/ilustre/users/sanger-dev/workspace/20190917/Single_findcirc2_8018_2388/Findcircbwa/output/rnasam.sam',
                'BSJ': '/mnt/ilustre/users/sanger-dev/workspace/20190927/Single_bsjreads_4702_6230/Bsjreads/output/zjx_circRNA_merge_signal_type_bsj.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


