# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class TrinityMergeAgent(Agent):
    """
    merge trinity partion result
    """
    def __init__(self, parent):
        super(TrinityMergeAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'partion_dir', 'format': 'denovo_rna_v2.common_dir'},
            {'type': 'outfile', 'name': 'trinity_fa', 'format': 'denovo_rna_v2.trinity_fasta'},
            {'default': 'denovo_rna_v2.common', 'type': 'outfile', 'name': 'gene2trans'},
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
        super(TrinityMergeAgent, self).end()


class TrinityMergeTool(Tool):
    """
    merge trinity partion result
    """
    def __init__(self, config):
        super(TrinityMergeTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.trinity_merge = software_dir + '/bioinfo/rna/scripts/trinity_merge.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_trinity_merge(self):
        cmd = '{} {} '.format(self.python_path, self.trinity_merge)
        cmd += '-{} {} '.format("partion_dir", self.option("partion_dir").prop['path'])
        cmd += '-{} {} '.format("trinity_fa", self.option("trinity_fa"))
        cmd += '-{} {} '.format("gene2trans", self.option("gene2trans"))
        cmd_name = 'trinity_merge'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32006501")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32006502")

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
        super(TrinityMergeTool, self).run()
        self.run_trinity_merge()
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
            "id": "TrinityMerge" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.trinity_merge",
            "instant": False,
            "options": dict(
                partion_dir=test_dir + "/" + "? infile name",
                trinity_fa="?",
                gene2trans="denovo_rna_v2.common",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
