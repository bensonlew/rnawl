# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd

class Mitarget2ceAgent(Agent):
    """
    Expression clustering analysis
    """
    def __init__(self, parent):
        super(Mitarget2ceAgent, self).__init__(parent)
        options = [
            {'name': 'mrna_target', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'lncrna_target', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'mirna_family', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'type', 'type': 'string', 'default': "G"},
            {'name': 'output', 'type': 'string', 'default': "ce.xls"},
            {'name': 'pvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("mrna_target").is_set:
            raise OptionError('mrna_target(small rna target file) should be set ')
        if not self.option("lncrna_target").is_set:
            raise OptionError('lncrna_target(lncrna target file) should be set ')

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
        super(Mitarget2ceAgent, self).end()


class Mitarget2ceTool(Tool):
    """
    Expression clustering analysis
    """
    def __init__(self, config):
        super(Mitarget2ceTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.ExpCorrsf_toolbox = self.config.PACKAGE_DIR + '/lnc_rna/target2ce.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)


    def run_target_ce(self):
        cmd = '{} {} '.format(self.python_path, self.ExpCorrsf_toolbox)
        cmd += '-{} {} '.format("mrna_target", self.option("mrna_target").prop['path'])
        cmd += '-{} {} '.format("lncrna_target", self.option("lncrna_target").prop['path'])
        cmd += '-{} {} '.format("mirna_family", self.option("mirna_family").prop['path'])
        cmd += '-{} {} '.format("gene_type", self.option("type"))

        '''
        if self.option("output") is None:
            self.option("output", self.work_dir + "corr.xls")
        else:
            if not os.path.exists(self.option("output")):
                os.mkdir(self.option("output"))
        '''
        cmd += '-{} {} '.format("output", self.option("output"))
        cmd_name = 'target_ce'
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
        super(Mitarget2ceTool, self).run()
        self.run_target_ce()
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
            "name": "lnc_rna.cerna.mitarget2ce",
            "instant": False,
            "options": dict(
                mrna_target = test_dir  + "/known_m.xls",
                lncrna_target = test_dir  + "/known_l.xls",
                mirna_family = test_dir  + "/known_miR_family.xls"

            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
