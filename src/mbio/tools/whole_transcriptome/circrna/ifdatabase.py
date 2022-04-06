# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd
import re

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class IfdatabaseAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(IfdatabaseAgent, self).__init__(parent)
        options = [
            {'name': 'details', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'circfasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            # {'name':'databasefasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'organism_name', 'type': 'string', 'default': None},
            {'name': 'detail_circbase', 'type': 'outfile', 'format': 'whole_transcriptome.common'},

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '15G'

    def end(self):
        super(IfdatabaseAgent, self).end()


class IfdatabaseTool(Tool):
    def __init__(self, config):
        super(IfdatabaseTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'basetomine': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/base2minetest.py')
        }
        self.file = {
            'detail_circbase': os.path.join(self.output_dir, 'detail.txt'),
            'database_path': os.path.join(self.config.SOFTWARE_DIR, 'database/Circrna/new/database_list')

        }

    def run(self):
        super(IfdatabaseTool, self).run()
        self.run_database()
        self.set_output()
        self.end()

    def run_database(self):
        database_path = read_database_path(self.file['database_path'], self.option('organism_name'))
        path = os.path.join(self.config.SOFTWARE_DIR, database_path)
        cmd = '{} {}'.format(self.program['python'], self.script['basetomine'])
        cmd += ' -i {} -d {} -m {} -o {}'.format(self.option('details').path, path, self.option('circfasta').path,
                                                 self.file['detail_circbase'])
        cmd_name = 'run_database'
        runcmd(self, cmd_name, cmd)

    # def run_database_id_new(self):
    #     circ = pd.read_table(self.file['detail_circbase'])
    #     database_id_new_list = list()
    #     for database_id in circ['database_id']:
    #         if database_id == '':
    #             database_id_new_list.append('')
    #         else:
    #             (chr, start, end) = re.split(r'[:-]', database_id)
    #             database_id_new = chr + ':' + str(int(start) + 1) + '-' + end
    #             database_id_new_list.append(database_id_new)
    #     circ['database_id'] = database_id_new_list
    #     circ.to_csv(self.file['detail_circbase'], index=False, header=True, sep='\t')

    def set_output(self):
        self.option('detail_circbase').set_path(self.file['detail_circbase'])


def read_database_path(database_list, organism_name):
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
            'id': 'ifdatabase_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.ifdatabase',
            'instant': False,
            'options': {
                'details': '/mnt/ilustre/users/sanger-dev/workspace/20191203/Circrna_tsg_36330/CircBrush/Circdetail/output/detail.txt',
                'circfasta': '/mnt/ilustre/users/sanger-dev/workspace/20191203/Circrna_tsg_36330/CircBrush/Getfasta/output/circrna.fasta',
                # 'databasefasta': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_getfasta_4419_1585/GetfastaCirc/output/circbase.fasta',
                'organism_name': 'Mus_musculus'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
