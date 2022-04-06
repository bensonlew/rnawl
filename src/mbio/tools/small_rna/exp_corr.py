# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import glob
from biocluster.tool import Tool
import pandas as pd
import unittest

class ExpCorrAgent(Agent):
    '''
    require express_matrix (file), scm (str), scd (str) and corr_method (str)
    '''
    def __init__(self, parent):
        super(ExpCorrAgent, self).__init__(parent)
        options = [
            {'name': 'express_matrix', 'type': 'infile', 'format': 'small_rna.express_matrix'},
            {'name': 'scm', 'type': 'string', 'default': 'complete'},
            {'name': 'scd', 'type': 'string', 'default': 'euclidean'},
            {'name': 'corr_method', 'type': 'string', 'default': 'pearson'},
            {'name': 'output', 'type': 'string', 'default': None},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option('express_matrix').prop['sample_number'] <= 1:
            raise OptionError('please input at least 2 samples')

    def set_resource(self):
        self._cpu = 2
        file_size = os.path.getsize(self.option('express_matrix').prop['path'])
        self._memory = '{}G'.format(int(float(file_size) / 1024 ** 3) + 5)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', 'exp_corr_tool_output_dir']
        ])
        super(ExpCorrAgent, self).end()

class ExpCorrTool(Tool):
    '''
    obtain exp_corr result, write data to sample_correlation.xls and expression_matrix.xls
    '''
    def __init__(self, config):
        super(ExpCorrTool, self).__init__(config)
        # old setting
        software_dir = self.config.SOFTWARE_DIR
        # self.python_path = 'miniconda2/bin/python'
        # self.cluster_toolbox = self.config.PACKAGE_DIR + '/denovo_rna_v2/cluster_toolbox.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        # self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        # self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        # self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        # self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

        # new setting, remove unnecessary definitions
        self.python_path = 'miniconda2/bin/python'
        self.cluster_toolbox = os.path.join(self.config.PACKAGE_DIR, 'small_rna/cluster_toolbox.py')
        self._PATH = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin')
        self._LD_LIBRARY_PATH = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/lib64/R/lib')
        self.set_environ(PATH=self._PATH, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_cluster_toolbox(self):
        exp_file = self.option('express_matrix').prop['path']
        t = pd.read_table(exp_file, index_col=0, header=0)
        if t.iloc[:,0].size > 200000:
            t1 = t[t.mean(1) > 0.5]
            exp_filter = os.path.join(self.work_dir, "exp_filter_matrix.xls")
            t1.to_csv(exp_filter, sep="\t")
        else:
            exp_filter = self.option('express_matrix').prop['path']
        cmd = '{} {} '.format(self.python_path, self.cluster_toolbox)
        cmd += '-exp {} '.format(exp_filter)
        # cmd += '-log_base 10 '
        cmd += '--ngc '
        cmd += '-{} {} '.format("scm", self.option("scm"))
        cmd += '-{} {} '.format("scd", self.option("scd"))
        cmd += '--corr '
        cmd += '-{} {} '.format("corr_method", self.option("corr_method"))
        if self.option('output') is None:
            self.option('output', self.work_dir)
        else:
            if not os.path.exists(self.option('output')):
                os.mkdir(self.option('output'))
        cmd += '-{} {} '.format('out', self.option('output'))
        cmd_name = 'exp_corr'
        command = self.add_command(cmd_name, cmd)
        self.logger.debug(cmd)
        self.logger.info('start command_run at tool')
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, return None, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}'.format(cmd_name))
        else:
            self.set_error('fail to run {}'.format(cmd_name))
        self.logger.info('finish command_run at tool')

    def set_output(self):
        self.logger.info('start set_output at tool')
        all_files = list()
        all_files.append(os.path.join(self.option('output'), 'sample_correlation.xls'))
        all_files.append(os.path.join(self.option('output'), 'sample.cluster_tree.txt'))
        all_files.append(os.path.join(self.option('output'), 'expression_matrix.xls'))
        if len(all_files) != 3:
            self.logger.debug('len of {} is {}'.format(all_files, len(all_files)))
            self.set_error('number of result files is incorrect')
        for source in all_files:
            basename = os.path.basename(source)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at tool')

    def run(self):
        super(ExpCorrTool, self).run()
        self.run_cluster_toolbox()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):

        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import random

        exp_matrix_list = ['/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_count.xls',
                           '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_norm.xls',
                           '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/novel_miR_count.xls',
                           '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/novel_miR_norm.xls']

        for exp_matrix in exp_matrix_list:

            data = {
                'id': 'exp_corr_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': 'tool',
                'name': 'small_rna.exp_corr',
                'instant': False,
                'options': {
                    'express_matrix': exp_matrix,
                    'scm': 'complete',
                    'scd': 'correlation',
                    'corr_method': 'pearson',
                    'output': None
                }
            }

            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

if __name__ == '__main__':
    unittest.main()
