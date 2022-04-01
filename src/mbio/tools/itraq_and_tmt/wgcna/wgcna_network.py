# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class WgcnaNetworkAgent(Agent):
    """
    wgcna network analysis 
    """
    def __init__(self, parent):
        super(WgcnaNetworkAgent, self).__init__(parent)
        options = [
            {'name': 'module', 'type': 'string'},
            {'name': 'threshold', 'type': 'string'},
            {'name': 'top', 'type': 'string', 'default': '30'},
            {'name': 'step3output', 'type': 'string'},
            {'name': 'step2output', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('30')

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
        super(WgcnaNetworkAgent, self).end()


class WgcnaNetworkTool(Tool):
    """
    wgcna network analysis 
    """
    def __init__(self, config):
        super(WgcnaNetworkTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.wgcna_network = self.config.PACKAGE_DIR + '/wgcna/wgcna_network.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1_gcc5.1/bin:$PATH"
        if "sanger-dev" in self.config.SOFTWARE_DIR:
            self.r_path = software_dir + "/program/R-3.3.3/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_wgcna_network(self):
        cmd = '{} {} '.format(self.python_path, self.wgcna_network)
        cmd += '-{} {} '.format("module", self.option("module"))
        cmd += '-{} {} '.format("threshold", self.option("threshold"))
        cmd += '-{} {} '.format("top", self.option("top"))
        cmd += '-{} {} '.format("step3output", self.option("step3output"))
        cmd += '-{} {} '.format("step2output", self.option("step2output"))
        cmd_name = 'wgcna_network'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708701")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708702")

    def set_output(self):
        all_files = glob.glob(self.work_dir + '/network.*')
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir,fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(WgcnaNetworkTool, self).run()
        self.run_wgcna_network()
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
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/rna/test_files'
        data = {
            "id": "WgcnaNetwork" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "rna.wgcna.wgcna_network",
            "instant": False,
            "options": dict(
                module="?",
                threshold="?",
                step3output="?",
                step2output="?",
                top="?"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
