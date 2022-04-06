# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd

class ExpcorrMitargetAgent(Agent):
    """
    Expression clustering analysis
    """
    def __init__(self, parent):
        super(ExpcorrMitargetAgent, self).__init__(parent)
        options = [
            {'name': 'exp', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'exp_target', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'rna_target', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'known_target', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'novol_target', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'output', 'type': 'string', 'default': "corr.xls"},
            {'name': 'pvalue_cutoff', 'type': 'float', 'default': 1.0},
            {'name': 'qvalue_cutoff', 'type': 'float', 'default': 1.0},
            {'name': 'corr_cutoff', 'type': 'float', 'default': 0},
            {'name': 'corr_way', 'type': 'string', 'default': "spearman"},
            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
            {'name': 'sig_type', 'type': 'int', 'default': 1},
        ]
        self.add_option(options)

    def check_options(self):

        if self.option("rna_target").is_set:
            pass
        elif self.option("known_target").is_set and self.option("novol_target").is_set:
            pass
        else:
            raise OptionError('rna_target or known/novol_target should be set ')

        if not self.option("exp").is_set:
            raise OptionError('exp(small rna expression file) should be set ')
        if not self.option("exp_target").is_set:
            raise OptionError('exp_target(target expression file) should be set ')

        return
        '''
        if self.option('exp').prop['gene_number'] >= 3001:
            raise OptionError('Please input at most 3000 genes', code = "33704701")
        '''

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
        super(ExpcorrMitargetAgent, self).end()


class ExpcorrMitargetTool(Tool):
    """
    Expression clustering analysis
    """
    def __init__(self, config):
        super(ExpcorrMitargetTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.ExpCorrsf_toolbox = self.config.PACKAGE_DIR + '/small_rna/exp_corr_mitarget.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def get_target_file(self, known_target, novol_target):
        with open("all_target.xls", "w") as all_w:
            with open(known_target, "r") as known_f:
                [all_w.write(line) for line in known_f]
            with open(novol_target, "r") as novol_f:
                novol_f.readline()
                [all_w.write(line) for line in novol_f]
        return "all_target.xls"

    def run_ExpCorrsf_toolbox(self):
        if self.option("rna_target").is_set:
            target = self.option("rna_target").prop['path']
        else:
            target = self.get_target_file(self.option("known_target").prop['path'], self.option("novol_target").prop['path'])
        cmd = '{} {} '.format(self.python_path, self.ExpCorrsf_toolbox)
        cmd += '-{} {} '.format("exp", self.option("exp").prop['path'])
        cmd += '-{} {} '.format("exp_target", self.option("exp_target").prop['path'])
        cmd += '-{} {} '.format("rna_target", self.option("rna_target").prop['path'])
        cmd += '-{} {} '.format("pvalue_cutoff", self.option("pvalue_cutoff"))
        cmd += '-{} {} '.format("qvalue_cutoff", self.option("qvalue_cutoff"))
        cmd += '-{} {} '.format("cor_cutoff", self.option("corr_cutoff"))
        cmd += '-{} {} '.format("corr_way", self.option("corr_way"))
        cmd += '-{} {} '.format("padjust_way", self.option("padjust_way"))

        '''
        if self.option("output") is None:
            self.option("output", self.work_dir + "corr.xls")
        else:
            if not os.path.exists(self.option("output")):
                os.mkdir(self.option("output"))
        '''
        cmd += '-{} {} '.format("output", self.option("output"))
        cmd_name = 'exp_corr_target'
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
                self.set_error("%s Failed. >>>%s"%(cmd_name, cmd) )
        else:
            self.set_error("%s Failed. >>>%s"%(cmd_name, cmd))

    def set_output(self):

        all_files = [self.option("output")]
        # all_files += glob.glob(self.option("output") + '/record.json')
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(ExpcorrMitargetTool, self).run()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA'
        data = {
            "id": "ExpCorrmitarget" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.expcorr_mitarget",
            "instant": False,
            "options": dict(
                exp=test_dir + "/data5" + "/known.norm.xls",
                exp_target = test_dir  + "/gene.fpkm.matrix.xls",
                rna_target = test_dir  + "/all_merged.xls",
                corr_way='pearson',
                output='corr.xls',
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
