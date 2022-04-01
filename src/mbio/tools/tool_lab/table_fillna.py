# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun,qinjincheng'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class TableFillnaAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(TableFillnaAgent, self).__init__(parent)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'sep', 'type': 'string'},
            {'name': 'nav', 'type': 'string'},
            {'name': 'rmp_type', 'type': 'bool'},
            {'name': 'rmp', 'type': 'float'},
            {'name': 'method', 'type': 'string'},
            {'name': 'transposed_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        super(TableFillnaAgent, self).end()


class TableFillnaTool(Tool):
    def __init__(self, config):
        super(TableFillnaTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
            'Rscript': os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/lib64/R/bin/Rscript')
        }
        self.script = {
            'table_kit': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/table_kit.py'),
            'fillna': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/fillna.r')
        }
        self.file = {
            'transposed_table': os.path.join(self.output_dir, '{}'.format(os.path.basename(self.option('table').path))),
        }

        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        # if "sanger-dev" in self.config.SOFTWARE_DIR:
        #     self.r_path = software_dir + "/program/R-3.5.1/bin:$PATH"
        # else:
        #     self.r_path = software_dir + "/program/R-3.5.1/bin:$PATH"
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        # self.r_path = software_dir + "/program/R-3.3.1_gcc5.1/bin:$PATH"
        # self._r_home = software_dir + "/program/R-3.5.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        # self._LD_LIBRARY_PATH = software_dir + "/program/R-3.5.1/lib64/R/library/Rcpp/libs:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        # self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
    def run(self):
        super(TableFillnaTool, self).run()

        self.run_fillna()
        self.set_output()
        self.end()

    def file_process(self):
        '''
        Remove blank lines
        Deal with '#' in header, as read.table() in R ignores lines that start with '#'
        '''
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[self.option('sep')]
        table_df = pd.read_table(self.option('table').prop['path'], sep=sep, index_col=None)
        table_df.dropna(how='all', inplace=True)
        if self.option('method') in ['missForest', 'kNN', 'PPCA', 'BPCA', 'nipals']:
            table_df.rename(columns={table_df.columns[0]: table_df.columns[0].replace('#', '')}, inplace=True)
        processed_file = os.path.join(self.work_dir, os.path.basename(self.option('table').prop['path']))
        table_df.to_csv(processed_file, sep=sep, index=None, header=True)
        return processed_file

    def run_fillna(self):
        input_table = self.file_process()   # added by zhangyitong on 20211202
        cmd = '{} {}'.format(self.program['python'], self.script['table_kit'])
        cmd += ' {}'.format('fillna')
        # cmd += ' --input {}'.format(self.option('table').path)
        cmd += ' --input {}'.format(input_table)
        cmd += ' --sep {}'.format(self.option('sep'))
        cmd += ' --nav {}'.format(self.option('nav'))
        if self.option('rmp_type') == True:
            cmd += ' --rmp {}'.format(self.option('rmp'))
        cmd += ' --method {}'.format(self.option('method'))
        cmd += ' --r_interpreter {}'.format(self.program['Rscript'])
        cmd += ' --r_script {}'.format(self.script['fillna'])
        cmd += ' --output {}'.format(self.file['transposed_table'])
        cmd_name = 'run_fillna'
        runcmd(self, cmd_name, cmd)


    def set_output(self):
        self.option('transposed_table').set_path(self.file['transposed_table'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'table_fillna_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.table_fillna',
            'instant': False,
            'options': {
                'table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/fillna_test.txt',
                'sep': 'tab',
                'nav': '0',
                'method': 'nipals',
                'rmp': '100',
                'rmp_type':True


            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


