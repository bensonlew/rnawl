# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd

class ExpCorrsfAgent(Agent):
    """
    Expression clustering analysis
    """
    def __init__(self, parent):
        super(ExpCorrsfAgent, self).__init__(parent)
        options = [
            {'name': 'exp', 'type': 'infile', 'format': 'denovo_rna_v2.express_matrix'},
            {'name': 'gt', 'type': 'string', 'default': None},
            {'name': 'anno', 'type': 'string', 'default': None},
            {'name': 'output', 'type': 'string', 'default': None},
            {'name': 'pvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'qvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'cor_cutoff', 'type': 'float', 'default': 0.8},
            {'name': 'corr_way', 'type': 'string', 'default': "spearmanr"},
            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
            {'name': 'sig_type', 'type': 'int', 'default': 1},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option('exp').prop['gene_number'] >= 3001:
            raise OptionError('Please input at most 3000 genes', code = "32007801")

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('20')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量相关分析结果目录"],
            ["./express_correlation_info.xls", "xls", "表达量相关性分析表"],
            ["./record.json", "json", "表达量相关性分析绘制网络图表"]
        ])
        super(ExpCorrsfAgent, self).end()


class ExpCorrsfTool(Tool):
    """
    Expression clustering analysis
    """
    def __init__(self, config):
        super(ExpCorrsfTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.ExpCorrsf_toolbox = self.config.PACKAGE_DIR + '/denovo_rna_v2/exp_corrsf.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_ExpCorrsf_toolbox(self):
        cmd = '{} {} '.format(self.python_path, self.ExpCorrsf_toolbox)
        cmd += '-{} {} '.format("exp", self.option("exp").prop['path'])
        cmd += '-{} {} '.format("g_or_t", self.option("gt"))
        cmd += '-{} {} '.format("anno", self.option("anno"))
        cmd += '-{} {} '.format("pvalue_cutoff", self.option("pvalue_cutoff"))
        cmd += '-{} {} '.format("qvalue_cutoff", self.option("qvalue_cutoff"))
        cmd += '-{} {} '.format("cor_cutoff", self.option("cor_cutoff"))
        cmd += '-{} {} '.format("corr_way", self.option("corr_way"))
        cmd += '-{} {} '.format("padjust_way", self.option("padjust_way"))
        cmd += '-{} {} '.format("sig_type", self.option("sig_type"))
        if self.option("output") is None:
            self.option("output", self.work_dir)
        else:
            if not os.path.exists(self.option("output")):
                os.mkdir(self.option("output"))
        cmd += '-{} {} '.format("output", self.option("output"))
        cmd_name = 'exp_corr'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32007801")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32007802")

    def set_output(self):

        all_files = glob.glob(self.option("output") + "/express_correlation_info.xls")
        all_files += glob.glob(self.option("output") + '/record.json')
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(ExpCorrsfTool, self).run()
        self.run_ExpCorrsf_toolbox()
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
            "id": "ExpCorrsf" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.exp_coranalysis",
            "instant": False,
            "options": dict(
                exp=test_dir + "/" + "transcript.tpm.matrix",
                scm="complete",
                scd="correlation",
                corr_method='pearson',
                output=None,
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
