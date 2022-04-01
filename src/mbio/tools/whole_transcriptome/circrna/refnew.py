# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class RefnewAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(RefnewAgent, self).__init__(parent)
        options = [
            {'name': 'circfasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'ref_no', 'type': 'outfile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new_no', 'type': 'outfile', 'format': 'whole_transcriptome.fasta'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '4G'

    def end(self):
        super(RefnewAgent, self).end()



class RefnewTool(Tool):
    def __init__(self, config):
        super(RefnewTool, self).__init__(config)
        self.file = {
            'ref_no': os.path.join(self.output_dir, 'ref.fasta'),
            'new_no': os.path.join(self.output_dir, 'new.fasta')


        }

    def run(self):
        super(RefnewTool, self).run()
        self.run_refnew()
        self.set_output()
        self.end()


    def run_refnew(self):
        fasta = open(self.option('circfasta').path, "r")
        fasta_all = fasta.read()
        handle_ref_no = open(self.file['new_no'], "w")
        handle_ref_no.write(fasta_all)
        handle_ref_no.close()
        handle_ref_new = open(self.file['ref_no'], "w")
        handle_ref_new.close()

    def set_output(self):
        self.option('ref_no').set_path(self.file['ref_no'])
        self.option('new_no').set_path(self.file['new_no'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'refnew_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.refnew',
            'instant': False,
            'options': {

                'circfasta': '/mnt/ilustre/users/sanger-dev/workspace/20191024/Single_circ_brush_4409_2962/CircBrush/output/circrna.fasta'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


