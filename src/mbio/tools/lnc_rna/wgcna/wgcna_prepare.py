# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class WgcnaPrepareAgent(Agent):
    """
    wgcna data pre-processing 
    """
    def __init__(self, parent):
        super(WgcnaPrepareAgent, self).__init__(parent)
        options = [
            {'name': 'exp', 'type': 'infile', 'format': 'rna.express_matrix'},
            {'name': 'me', 'type': 'float', 'default': 0},
            {'name': 'cv', 'type': 'float', 'default': 0},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 3
        self._memory = "{}G".format(30)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(WgcnaPrepareAgent, self).end()


class WgcnaPrepareTool(Tool):
    """
    wgcna data pre-processing 
    """
    def __init__(self, config):
        super(WgcnaPrepareTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.wgcna_prepare = self.config.PACKAGE_DIR + '/wgcna/wgcna_prepare.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        if "sanger-dev" in self.config.SOFTWARE_DIR:
            self.r_path = software_dir + "/program/R-3.3.3/bin:$PATH"
        else:
            self.r_path = software_dir + "/program/R-3.3.1_gcc5.1/bin:$PATH"
        # self.r_path = software_dir + "/program/R-3.3.1_gcc5.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_wgcna_prepare(self):
        cmd = '{} {} '.format(self.python_path, self.wgcna_prepare)
        cmd += '-{} {} '.format("exp", self.option("exp").prop['path'])
        cmd += '-{} {} '.format("me", self.option("me"))
        cmd += '-{} {} '.format("cv", self.option("cv"))
        cmd_name = 'wgcna_prepare'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708801")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708802")

    def set_output(self):
        f1 = glob.glob(self.work_dir + '/scale_free_analysis.xls')
        f2 = glob.glob(self.work_dir + '/sample.cluster.dendrogram.*')
        f3 = glob.glob(self.work_dir + '/powerEstimate_*')
        f4 = glob.glob(self.work_dir + '/ignored_gene.list')
        f5 = glob.glob(self.work_dir + '/exp_matrix_after_filtering.txt')
        all_files = f1 + f2 + f3 + f4 + f5
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(WgcnaPrepareTool, self).run()
        self.run_wgcna_prepare()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "WgcnaPrepare" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "rna.wgcna.wgcna_prepare",
            "instant": False,
            "options": dict(
                exp="/mnt/ilustre/users/sanger-dev/sg-users/deqing/wgcna_test2/exp_matrix.txt",
                me="0.1",
                cv="0.1",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
