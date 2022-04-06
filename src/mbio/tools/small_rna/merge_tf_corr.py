# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'liubinxu'


class MergeTfCorrAgent(Agent):
    """
    merge tf binding prediction and corr of mirna_target
    """
    def __init__(self, parent):
        super(MergeTfCorrAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'tfbinding', 'format': 'small_rna.common'},
            {'type': 'infile', 'name': 'corr', 'format': 'small_rna.common'},
            {'name': 'nodes', 'type': 'string', 'default': 'nodes.xls'},
            {'name': 'edges', 'type': 'string', 'default': 'edges.xls'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('10')

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
        super(MergeTfCorrAgent, self).end()


class MergeTfCorrTool(Tool):
    """
    merge tf binding prediction and corr of mirna_target
    """
    def __init__(self, config):
        super(MergeTfCorrTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.python_path = 'miniconda2/bin/python'
        self.net_merge = self.config.PACKAGE_DIR + '/small_rna/net_merge.py'

    def run_net_merge(self):
        if self.option("tfbinding").is_set:
            tf_file = self.option("tfbinding").prop['path']
        else:
            tf_file = None
        cmd = '{} {} '.format(self.python_path, self.net_merge)
        cmd += ' {} '.format(tf_file)
        cmd += ' {} '.format(self.option("corr").prop['path'])
        cmd += ' {} '.format(self.option("nodes"))
        cmd += ' {} '.format(self.option("edges"))
        cmd_name = 'merge_tf_corr'
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
        pass
        # '''Example:
        all_files = [self.option("nodes"), self.option("edges")]
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)


    def run(self):
        super(MergeTfCorrTool, self).run()
        self.run_net_merge()
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
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA'
        
        data = {
            "id": "MergeTfCorr" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.merge_tf_corr",
            "instant": False,
            "options": dict(
                tfbinding=test_dir + "/binding.xls",
                corr=test_dir + "/corr.xls",
                nodes="nodes.xls",
                edges="net.xls",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
