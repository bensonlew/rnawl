# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd

class ExpcorrCernaAgent(Agent):
    """
    Expression clustering analysis
    """
    def __init__(self, parent):
        super(ExpcorrCernaAgent, self).__init__(parent)
        options = [
            {'name': 'exp_mirna', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'exp_cerna', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'corr_mirna2mrna', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'corr_mirna2lncrna', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ce_pair', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'output', 'type': 'string', 'default': "corr.xls"},
            {'name': 'pvalue_cutoff', 'type': 'float', 'default': 1},
            {'name': 'qvalue_cutoff', 'type': 'float', 'default': 1},
            {'name': 'corr_cutoff', 'type': 'float', 'default': 0.0},
            {'name': 'corr_way', 'type': 'string', 'default': "spearman"},
            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
            {'name': 'type', 'type': 'string', 'default': 'G'},
            {'name': 'sig_type', 'type': 'int', 'default': 1},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("ce_pair").is_set and self.option("exp_cerna").is_set:
            pass
        else:
            raise OptionError('ce_pair or exp_cerna should be set ')
        return


    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('20')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量相关分析结果目录"],
            ["./express_correlation_info.xls", "xls", "表达量相关性分析表"],
        ])
        super(ExpcorrCernaAgent, self).end()


class ExpcorrCernaTool(Tool):
    """
    Expression clustering analysis
    """
    def __init__(self, config):
        super(ExpcorrCernaTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.ExpCorrsf_toolbox = self.config.PACKAGE_DIR + '/lnc_rna/exp_corr_ce.py'
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
        '''
        if self.option("rna_target").is_set:
            target = self.option("rna_target").prop['path']
        else:
            target = self.get_target_file(self.option("known_target").prop['path'], self.option("novol_target").prop['path'])
        '''


        cmd = '{} {} '.format(self.python_path, self.ExpCorrsf_toolbox)
        if self.option("exp_mirna").is_set:
            cmd += '-{} {} '.format("exp_mirna", self.option("exp_mirna").prop['path'])
        else:
            pass
        cmd += '-{} {} '.format("exp_cerna", self.option("exp_cerna").prop['path'])
        cmd += '-{} {} '.format("ce_pair", self.option("ce_pair").prop['path'])
        cmd += '-{} {} '.format("pvalue_cutoff", self.option("pvalue_cutoff"))
        cmd += '-{} {} '.format("qvalue_cutoff", self.option("qvalue_cutoff"))
        cmd += '-{} {} '.format("cor_cutoff", self.option("corr_cutoff"))
        cmd += '-{} {} '.format("corr_way", self.option("corr_way"))
        cmd += '-{} {} '.format("padjust_way", self.option("padjust_way"))
        cmd += '-{} {} '.format("g_or_t", self.option("type"))
        cmd += '-{} {} '.format("output", self.option("output"))

        corr_m2ce_list = []
        if self.option("corr_mirna2mrna").is_set:
            corr_m2ce_list.append(self.option("corr_mirna2mrna").prop['path'])
        if self.option("corr_mirna2lncrna").is_set:
            corr_m2ce_list.append(self.option("corr_mirna2lncrna").prop['path'])

        cmd += '-{} {} '.format("corr_mirna2cerna", ",".join(corr_m2ce_list))
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
        super(ExpcorrCernaTool, self).run()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna'
        data = {
            "id": "ExpCorrmitarget" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "lnc_rna.cerna.expcorr_cerna",
            "instant": False,
            "options": dict(
                exp_mirna=test_dir + "/mirna.exp.xls",
                exp_cerna = test_dir  + "/ce.exp.xls",
                corr_mirna2mrna = test_dir  + "/known_m.corr.xls",
                corr_mirna2lncrna = test_dir  + "/known_l.corr.xls",
                ce_pair = test_dir + "/ce.xls",
                corr_way='pearson',
                output='corr.xls',
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
