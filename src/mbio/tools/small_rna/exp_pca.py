# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest

class ExpPcaAgent(Agent):
    """
    Expression clustering analysis
    """
    def __init__(self, parent):
        super(ExpPcaAgent, self).__init__(parent)
        options = [
            {'name': 'express_matrix', 'type': 'infile', 'format': 'small_rna.express_matrix'},
            {'name': 'output', 'type': 'string', 'default': None},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option('express_matrix').prop['sample_number'] <= 1:
            raise OptionError('Please input at least 2 samples')

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('5')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "exp_pca_tool_output_dir"]
            ])
        super(ExpPcaAgent, self).end()


class ExpPcaTool(Tool):
    """
    Expression clustering analysis
    """
    def __init__(self, config):
        super(ExpPcaTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        # self.cluster_toolbox = software_dir + '/bioinfo/rna/scripts/cluster_toolbox.py'
        self.cluster_toolbox = self.config.PACKAGE_DIR + '/small_rna/cluster_toolbox.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_cluster_toolbox(self):
        cmd = '{} {} '.format(self.python_path, self.cluster_toolbox)
        cmd += '-exp {} '.format(self.option('express_matrix').prop['path'])
        cmd += '-log_base 10 '
        cmd += '--ngc '
        cmd += '--nsc '
        cmd += '--pca '
        if self.option("output") is None:
            self.option("output", self.work_dir)
        else:
            if not os.path.exists(self.option("output")):
                os.mkdir(self.option("output"))
        cmd += '-{} {} '.format("out", self.option("output"))
        cmd_name = 'exp_pca'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32003002")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32003003")

    def set_output(self):
        all_files = glob.glob(self.option("output") + '/PCA.xls')
        all_files += glob.glob(self.option("output") + '/Explained_variance_ratio.xls')
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(ExpPcaTool, self).run()
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
                'id': 'exp_pca_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': 'tool',
                'name': 'small_rna.exp_pca',
                'instant': False,
                'options': {
                    'express_matrix': exp_matrix,
                    'output': None
                }
            }

            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

if __name__ == '__main__':
    unittest.main()
