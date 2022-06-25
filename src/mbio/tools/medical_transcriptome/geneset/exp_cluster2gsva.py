# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import pandas as pd
__author__ = 'gdq'


class ExpCluster2gsvaAgent(Agent):
    """
    Expression clustering analysis
    """
    def __init__(self, parent):
        super(ExpCluster2gsvaAgent, self).__init__(parent)
        options = [
            {'name': 'exp', 'type': 'string'},
            {'name': 'group', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'use_group', 'type': 'string', 'default': 'no'},
            {'name': 'sct', 'type': 'string', 'default': 'hierarchy'},
            {'name': 'gct', 'type': 'string', 'default': 'hierarchy'},
            {'name': 'scm', 'type': 'string', 'default': 'complete'},
            {'name': 'gcm', 'type': 'string', 'default': 'average'},
            {'name': 'scd', 'type': 'string', 'default': 'correlation'},
            {'name': 'gcd', 'type': 'string', 'default': 'euclidean'},
            {'name': 'output', 'type': 'string', 'default': None},
            {'name': 'static', 'type': 'bool', 'default': True},
            {'name': 'n_clusters', 'type': 'string'}
        ]
        self.add_option(options)

    def check_options(self):
        if self.option('sct') == 'hclust':
            self.option('sct', 'hierarchy')
        if self.option('gct') == 'hclust':
            self.option('gct', 'hierarchy')
        self._memory_increase_step = 30

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
        super(ExpCluster2gsvaAgent, self).end()


class ExpCluster2gsvaTool(Tool):
    """
    Expression clustering analysis
    """
    def __init__(self, config):
        super(ExpCluster2gsvaTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.cluster_toolbox = self.config.PACKAGE_DIR + '/medical_transcriptome/geneset/cluster_toolbox2gsva.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/bioinfo/miniconda2/bin:$PATH"
        self._r_home = software_dir + "/bioinfo/miniconda2/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/bioinfo/miniconda2/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.R = software_dir + "/bioinfo/miniconda2/bin/Rscript"
        self.Rscrip = self.config.PACKAGE_DIR + '/denovo_rna_v2/pheatmap.r'
    def run_cluster_toolbox(self):
        df = pd.read_table(self.option('exp'))
        # if len(df) > 3000:
        #     self.option('static', True)
        # else:
        #     self.option('static', False)
        self.option('static', False)
        cmd = '{} {} '.format(self.python_path, self.cluster_toolbox)
        cmd += '-{} {} '.format("exp", self.option("exp"))
        if self.option('use_group') != 'no':
            cmd += '-{} {} '.format("group", self.option("group").prop['path'])
        # cmd += '-{} {} '.format("n_clusters", None)
        if self.option('sct') == 'no':
            cmd += '--nsc '
        if self.option('gct') == 'no':
            cmd += '--ngc '
        cmd += '-{} {} '.format("sct", self.option("sct"))
        cmd += '-{} {} '.format("gct", self.option("gct"))
        cmd += '-{} {} '.format("scm", self.option("scm"))
        cmd += '-{} {} '.format("gcm", self.option("gcm"))
        cmd += '-{} {} '.format("scd", self.option("scd"))
        cmd += '-{} {} '.format("gcd", self.option("gcd"))
        if self.option('static'):
            cmd += '-{} {} '.format('static', self.option('static'))
            cmd += '-{} {} '.format('R', self.R)
            cmd += '-{} {} '.format('Rscrip', self.Rscrip)
            cmd += '-{} {} '.format('group4heatmap', self.option('group').prop['path'])
        if self.option("output") is None:
            self.option("output", self.work_dir)
        else:
            if not os.path.exists(self.option("output")):
                os.mkdir(self.option("output"))
        cmd += '-{} {} '.format("out", self.option("output"))
        cmd_name = 'exp_cluster'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33704402")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33704403")

    def set_output(self):

        all_files = glob.glob(self.option("output") + '/sample*cluster*tree.txt')
        all_files += glob.glob(self.option("output") + '/seq*cluster*.txt')
        all_files += glob.glob(self.option("output") + '/seq*cluster*.xls')
        all_files += glob.glob(self.option("output") + '/expression_matrix.xls')
        all_files += glob.glob(self.option('output') + '/heatmap.pdf')
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(ExpCluster2gsvaTool, self).run()
        self.run_cluster_toolbox()
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
            "id": "ExpCluster" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.exp_cluster",
            "instant": False,
            "options": dict(
                exp="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/exp_matrix",
                group="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/group",
                sct="hierarchy",
                gct="hierarchy",
                scm="complete",
                gcm="complete",
                scd="euclidean",
                gcd="euclidean",
                # sct="no",
                # gct="no",
                # static=True
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
