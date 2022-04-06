# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class SplitFastaAgent(Agent):
    """
    split_fasta description
    """
    def __init__(self, parent):
        super(SplitFastaAgent, self).__init__(parent)
        options = [
            {'name': 'f', 'type': 'infile', 'default': '/mnt/ilustre/users/isanger/sg-users/deqing/TestFiles/TF/all_pep.fa', 'format': 'sequence.fasta'},
            {'name': 'size', 'type': 'int', 'default': '1000', 'format': 'None'},
            {'name': 'prefix', 'type': 'string', 'default': 'split', 'format': 'None'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('15')

    def end(self):
        super(SplitFastaAgent, self).end()


class SplitFastaTool(Tool):
    """
    split_fasta description
    """
    def __init__(self, config):
        super(SplitFastaTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.split_fasta = self.config.PACKAGE_DIR + '/transcription_factor/split_fasta.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_split_fasta(self):
        cmd = '{} {} '.format(self.python_path, self.split_fasta)
        cmd += '-{} {} '.format("f", self.option("f").prop['path'])
        cmd += '-{} {} '.format("size", self.option("size"))
        cmd += '-{} {} '.format("prefix", self.option("prefix"))
        cmd_name = 'split_fasta'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708101")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708102")

    def set_output(self):
        target_files = glob.glob(self.work_dir + '/{}*.fa'.format(self.option('prefix')))
        for each in target_files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(SplitFastaTool, self).run()
        self.run_split_fasta()
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
            "id": "SplitFasta" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "rna.split_fasta",
            "instant": False,
            "options": dict(
                f="/mnt/ilustre/users/isanger/sg-users/deqing/TestFiles/TF/all_pep.fa",
                size="1000",
                prefix="split",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
