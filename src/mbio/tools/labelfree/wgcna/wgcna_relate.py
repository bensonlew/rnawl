# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class WgcnaRelateAgent(Agent):
    """
    wgcna relate analysis 
    """
    def __init__(self, parent):
        super(WgcnaRelateAgent, self).__init__(parent)
        options = [
            {'name': 'datExpr', "type": "infile", "format": "labelfree.common"},
            {'name': 'MEs', 'type': 'string'},
            {'name': 'traits', 'type': 'string'},
            {'name': 'corType', 'type': 'string', 'default': 'pearson'},
            {'name': 'nThreads', 'type': 'int', 'default': 16},
            {'name': 'block_Rdata', 'type': 'string', 'default': None},
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
        super(WgcnaRelateAgent, self).end()


class WgcnaRelateTool(Tool):
    """
    wgcna relate analysis 
    """
    def __init__(self, config):
        super(WgcnaRelateTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.wgcna_relate = self.config.PACKAGE_DIR + '/wgcna/wgcna_relate.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1_gcc5.1/bin:$PATH"
        if "sanger-dev" in self.config.SOFTWARE_DIR:
            self.r_path = software_dir + "/program/R-3.3.3/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_wgcna_relate(self):
        cmd = '{} {} '.format(self.python_path, self.wgcna_relate)
        cmd += '-{} {} '.format("datExpr", self.option("datExpr").prop['path'])
        cmd += '-{} {} '.format("MEs", self.option("MEs"))
        cmd += '-{} {} '.format("traits", self.option("traits"))
        cmd += '-{} {} '.format("corType", self.option("corType"))
        cmd += '-{} {} '.format("nThreads", self.option("nThreads"))
        if self.option("block_Rdata"):
            cmd += '-{} {} '.format("block_Rdata", self.option("block_Rdata"))
        cmd_name = 'wgcna_relate'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708901")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708902")

    def set_output(self):
        all_files = glob.glob(self.work_dir + '/module_trait.correlation.xls')
        all_files += glob.glob(self.work_dir + '/module_trait.correlation_pvalues.xls')
        all_files += glob.glob(self.work_dir + '/gene_trait.correlation.xls')
        all_files += glob.glob(self.work_dir + '/*gene_trait.correlation.pdf')
        all_files += glob.glob(self.work_dir + '/*gene_trait.correlation.png')
        for each in all_files:
            if 'gene' in each:
                each_ = each.replace('gene', 'protein')
                os.rename(each, each_)
                each = each_
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(WgcnaRelateTool, self).run()
        self.run_wgcna_relate()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/deqing/mbio/tools/rna/wgcna/test_files'
        data = {
            "id": "WgcnaRelate" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "rna.wgcna.wgcna_relate",
            "instant": False,
            "options": dict(
                datExpr=test_dir + "/" + "exp_matrix_after_filtering.txt",
                MEs=test_dir + "/" + "eigengenes.txt",
                traits=test_dir + "/" + "traits.xls",
                corType="pearson",
                nThreads="16",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
