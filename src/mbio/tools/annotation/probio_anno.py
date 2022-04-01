# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import unittest


class ProbioAnnoAgent(Agent):
    """
    宏基因组根据NR的gene_anno注释结果和益生菌数据库比对，得到益生菌注释结果表
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(ProbioAnnoAgent, self).__init__(parent)
        options = [
            {"name": "nr_gene_anno", "type": "infile", "format": "sequence.profile_table", "required": True},
            # nr的基因注释结果表
            {"name": "probio_anno", "type": "outfile", 'format': "sequence.profile_table"}  # 注释详细结果表
        ]
        self.add_option(options)

    def check_options(self):
        '''
        '''
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        self.option('probio_anno', os.path.join(self.output_dir, "gene_probio_anno.xls"))
        super(ProbioAnnoAgent, self).end()


class ProbioAnnoTool(Tool):
    def __init__(self, config):
        super(ProbioAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/probio.py'
        self.ref_file = self.config.SOFTWARE_DIR + "/database/Probio/Probio_for_all.txt"

    def run(self):
        """
        运行
        :return:
        """
        super(ProbioAnnoTool, self).run()
        self.run_probio_anno()
        self.end()

    def run_probio_anno(self):
        out_file = self.output_dir + "/gene_probio_anno.xls"
        gene_anno_nr = self.option("nr_gene_anno").prop["path"]
        cmd = '{} {} -i {} -ref {} -o {}'.format(self.python_path, self.python_script, gene_anno_nr, self.ref_file, out_file)
        self.logger.info(cmd)
        #self.wait(command)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('run probio_anno succeed')
        except subprocess.CalledProcessError:
            self.set_error("probio_anno failed", code="31205101")
            self.set_error('probio_anno failed', code="31205102")



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "probio_anno",
            "type": "tool",
            "name": "annotation.probio_anno",
            "instant": True,
            "options": dict(
                nr_gene_anno="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/metag_v2/nr/gene_nr_anno.xls",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
