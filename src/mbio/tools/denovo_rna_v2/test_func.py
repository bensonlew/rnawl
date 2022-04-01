# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class TestFuncAgent(Agent):
    """
    This is tool description
    """
    def __init__(self, parent):
        super(TestFuncAgent, self).__init__(parent)
        options = [
            {'name': 'arg1', 'type': 'float', 'default': 0.05},
            {'name': 'arg2', 'type': 'infile', 'default': 3, 'format': 'denovo_rna_v2.express_matrix'},
            {'name': 'arg3', 'type': 'string', 'default': '4'},
            {'name': 'arg4', 'type': 'int', 'default': 0},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('?')

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
        super(TestFuncAgent, self).end()


class TestFuncTool(Tool):
    """
    This is tool description
    """
    def __init__(self, config):
        super(TestFuncTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.called_script_name = software_dir + '/bioinfo/rna/scripts/called_script_name.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_called_script_name(self):
        cmd = '{} {} '.format(self.python_path, self.called_script_name)
        cmd += '-{} {} '.format("arg1", self.option("arg1"))
        cmd += '-{} {} '.format("arg2", self.option("arg2").prop['path'])
        cmd += '-{} {} '.format("arg3", self.option("arg3"))
        cmd += '-{} {} '.format("arg4", self.option("arg4"))
        cmd_name = 'TestFunc'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32005801")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32005802")

    def set_output(self):
        pass
        '''Example:
        diff_files = glob.glob(self.option("output") + '/*_vs_*.xls')
        diff_list = glob.glob(self.option("output") + '/*.DE.list')
        diff_summary = glob.glob(self.option("output") + '/*summary.xls')
        all_files = diff_files + diff_list + diff_summary
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
        '''
    def run(self):
        super(TestFuncTool, self).run()
        self.run_TestFunc()
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
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'
        data = {
            "id": "TestFunc" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.test_func",
            "instant": True,
            "options": dict(
                arg1="0.05",
                arg2="3",
                arg3="4",
                arg4="0",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
