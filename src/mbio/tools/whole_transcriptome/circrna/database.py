# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class DatabaseAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(DatabaseAgent, self).__init__(parent)
        options = [
            {'name': 'details', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            # {'name':'database', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'organism_name', 'type': 'string', 'default': None},
            {'name': 'detail_circbase', 'type': 'outfile', 'format':'whole_transcriptome.common'}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '4G'

    def end(self):
        super(DatabaseAgent, self).end()



class DatabaseTool(Tool):
    def __init__(self, config):
        super(DatabaseTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'database': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/database.py')
        }
        self.file = {
            'detail_circbase': os.path.join(self.output_dir, 'detail.txt'),
            'database_path': os.path.join(self.config.SOFTWARE_DIR, 'database/Circrna/database_list')

        }

    def run(self):
        super(DatabaseTool, self).run()
        self.run_database()
        self.set_output()
        self.end()


    def run_database(self):
        database_path = read_database_path(self.file['database_path'],self.option('organism_name'))
        cmd = '{} {}'.format(self.program['python'], self.script['database'])
        cmd += ' -i {} -d {} -r {} -o {}'.format(self.option('details').path, database_path, self.option('organism_name'),self.file['detail_circbase'])
        cmd_name = 'run_database'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('detail_circbase').set_path(self.file['detail_circbase'])

def read_database_path(database_list,organism_name):
    for line in open(database_list):
        eles = line.rstrip().split(' ')
        database_path = eles[1]
        organism = eles[0]
        if organism == organism_name:
            return database_path

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'database_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.database',
            'instant': False,
            'options': {
                'details':'/mnt/ilustre/users/sanger-dev/workspace/20191016/Single_circ_brush_6168_2072/CircBrush/output/detail.txt',
                # 'database':'/mnt/ilustre/users/sanger-dev/app/database/Circrna/Homo_sapiens_hg19_to_GRCh38.txt',
                'organism_name':'Homo_sapiens'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


