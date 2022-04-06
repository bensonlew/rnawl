# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'liubinxu'


class Ncbi2cogAgent(Agent):
    """
    convert Ncbi to cog
    """
    def __init__(self, parent):
        super(Ncbi2cogAgent, self).__init__(parent)
        options = [
            {'default': '', 'type': 'string', 'name': 'blasttable'},
            {'default': 'cog.xls', 'type': 'string', 'name': 'cogtable'},
            {"name": "cog_version", "type": "string", "default": "2020"},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option('blasttable') == "":
            raise OptionError("必须提供BLAST结果文件", code = "35000501")
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('5')

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
        super(Ncbi2cogAgent, self).end()


class Ncbi2cogTool(Tool):
    """
    convert Ncbi to cog
    """
    def __init__(self, config):
        super(Ncbi2cogTool, self).__init__(config)
        package_dir = self.config.PACKAGE_DIR
        self.python = "miniconda2/bin/python"
        self.mg_cog_mongo = package_dir + '/prok_rna/cog_ncbi.py'


    def run_mg_cog_mongo(self):
        cmd = '{} {} '.format(self.python, self.mg_cog_mongo)
        cmd += '-{} {} '.format("i", self.option("blasttable"))
        cmd += '-{} {} '.format("o", self.option("cogtable"))
        cmd += '-{} {} '.format("v", self.option("cog_version"))
        cmd_name = 'ncbi2cog'
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
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def set_output(self):
        basename = os.path.basename(self.option("cogtable"))
        link = os.path.join(self.output_dir, basename)
        if os.path.exists(link):
            os.remove(link)
        os.link(self.option("cogtable"), link)

    def run(self):
        super(Ncbi2cogTool, self).run()
        self.run_mg_cog_mongo()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_prok_rna'
        data = {
            "id": "eggnog2cog" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna.annotation.eggnog2cog",
            "instant": False,
            "options": dict(
                blasttable= test_dir + "/annotation/blast.xls",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
