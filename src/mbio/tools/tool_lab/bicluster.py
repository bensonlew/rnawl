# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class BiclusterAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(BiclusterAgent, self).__init__(parent)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "preprocess", "type": "string", 'default': 'no'},
            {"name": "format", "type": "string", 'default': "read_counts"},
            {"name": "missing", "type": "string", 'default': "geneMedian"},
            {"name": "ngenes", "type": "int", 'default': 2000},
            {"name": "method", "type": "string", 'default': 'BCCC'},
        ]
        self.add_option(options)
        self.step.add_steps('bicluster')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.bicluster.start()
        self.step.update()

    def step_end(self):
        self.step.bicluster.finish()
        self.step.update()

    def check_options(self):
        if not self.option('infile').is_set:
            raise OptionError('必须输入表达量文件。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(BiclusterAgent, self).end()


class BiclusterTool(Tool):
    def __init__(self, config):
        super(BiclusterTool, self).__init__(config)
        # self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        # self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.rscript = 'bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/Rscript'
        self.bicluster_R = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bicluster.r')

    def do_bicluster(self):
        """
        """
        cmd = "{} {} ".format(self.rscript, self.bicluster_R)
        cmd += '--infile {} '.format(self.option('infile').prop['path'])
        cmd += '--outdir {} '.format(self.output_dir)
        if self.option('preprocess') == "yes":
            cmd += '--preprocess '
        if self.option('format') == "read_counts":
            cmd += '--format 1 '
        if self.option('format') == "normalized_expression":
            cmd += '--format 2 '
        cmd += '--missing {} '.format(self.option('missing'))
        #
        cmd += '--ngenes {} '.format(self.option('ngenes'))
        cmd += "--method '{}()' ".format(self.option('method'))
        cmd += '--heatmap '
        self.logger.info(cmd)
        self.logger.info("开始进行双聚类bicluster")
        command = self.add_command("bicluster", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("双聚类bicluster完成")
        else:
            self.set_error("双聚类bicluster出错！")

    def run(self):
        super(BiclusterTool, self).run()
        self.do_bicluster()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "bicluster_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.bicluster",
            "instant": False,
            "options": dict(
                infile='/mnt/ilustre/users/sanger-dev/sg-users/xuxi/bicluster_test_dir/counts_cmp_use_small_good.txt',
                preprocess='yes',
                format="read_counts",
                missing="geneMedian",
                ngenes=2000,
                method="BCCC"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)