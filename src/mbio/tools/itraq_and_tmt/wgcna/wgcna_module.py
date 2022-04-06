# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class WgcnaModuleAgent(Agent):
    """
    wgcna module identification
    """
    def __init__(self, parent):
        super(WgcnaModuleAgent, self).__init__(parent)
        options = [
            {'name': 'datExpr', 'type': 'infile', 'format': 'rna.express_matrix'},
            {'name': 'mergeCutHeight', 'type': 'float', 'default': 0.25},
            {'name': 'corType', 'type': 'string', 'default': 'pearson'},
            {'name': 'maxBlockSize', 'type': 'int', 'default': 20000},
            {'name': 'power', 'type': 'int', 'default': 6},
            {'name': 'networkType', 'type': 'string', 'default': 'signed'},
            {'name': 'minModuleSize', 'type': 'int', 'default': 0},
            {'name': 'minKMEtoStay', 'type': 'float', 'default': 0.3},
            {'name': 'nThreads', 'type': 'int', 'default': 16},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
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
        super(WgcnaModuleAgent, self).end()


class WgcnaModuleTool(Tool):
    """
    wgcna module identification
    """
    def __init__(self, config):
        super(WgcnaModuleTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.wgcna_module = self.config.PACKAGE_DIR + '/wgcna/wgcna_module.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1_gcc5.1/bin:$PATH"
        if "sanger-dev" in self.config.SOFTWARE_DIR:
            self.r_path = software_dir + "/program/R-3.3.3/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_wgcna_module(self):
        cmd = '{} {} '.format(self.python_path, self.wgcna_module)
        cmd += '-{} {} '.format("datExpr", self.option("datExpr").prop['path'])
        cmd += '-{} {} '.format("mergeCutHeight", self.option("mergeCutHeight"))
        cmd += '-{} {} '.format("corType", self.option("corType"))
        cmd += '-{} {} '.format("maxBlockSize", self.option("maxBlockSize"))
        cmd += '-{} {} '.format("power", self.option("power"))
        cmd += '-{} {} '.format("networkType", self.option("networkType"))
        cmd += '-{} {} '.format("minModuleSize", self.option("minModuleSize"))
        cmd += '-{} {} '.format("minKMEtoStay", self.option("minKMEtoStay"))
        cmd += '-{} {} '.format("nThreads", self.option("nThreads"))
        cmd_name = 'wgcna_module'
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
        f1 = glob.glob(self.work_dir + '/module_corr.matrix.xls')
        f2 = glob.glob(self.work_dir + '/membership.xls')
        f3 = glob.glob(self.work_dir + '/module_size.stat.xls')
        f4 = glob.glob(self.work_dir + '/eigen*s.txt')
        f5 = glob.glob(self.work_dir + '/TOM-*.RData')
        f6 = glob.glob(self.work_dir + '/block_*_dendrogram.*pdf')
        f7 = glob.glob(self.work_dir + '/blockwiseModules_result.RData')
        f8 = glob.glob(self.work_dir + '/../seq_id2*_name.txt')
        f9 = glob.glob(self.work_dir + '/block_*_dendrogram.*png')
        all_files = f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9
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
        super(WgcnaModuleTool, self).run()
        self.run_wgcna_module()
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
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/rna/wgcna/test_files'
        data = {
            "id": "WgcnaModule" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "rna.wgcna.wgcna_module",
            "instant": False,
            "options": dict(
                datExpr=test_dir + "/" + "exp_matrix.txt",
                mergeCutHeight="0.25",
                corType="pearson",
                maxBlockSize="20000",
                power="6",
                networkType="signed",
                minModuleSize="0",
                minKMEtoStay="0.3",
                nThreads="16",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
