# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import unittest


class ProbioAbunAgent(Agent):
    """
    宏基因组益生菌
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(ProbioAbunAgent, self).__init__(parent)
        options = [
            {"name": "probio_anno", "type": "infile", "format": "sequence.profile_table", "required": True},
            # probio的基因注释结果表
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table", "required": True},
            # 基因丰度表
            {"name": "probio_abu", "type": "outfile", 'format': "sequence.profile_table"}  # 注释详细结果表
        ]
        self.add_option(options)

    def check_options(self):
        '''
        if not self.option("nr_gene_anno").is_set:
            raise OptionError("必须设置输入文件", code="")
        '''
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        self.option('probio_abu', os.path.join(self.output_dir, "probio_abun.xls"))
        super(ProbioAbunAgent, self).end()


class ProbioAbunTool(Tool):
    def __init__(self, config):
        super(ProbioAbunTool, self).__init__(config)
        self._version = "1.0"
        #self.soft = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl' 
        self.soft = '/miniconda2/bin/perl'
        self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/probio_abu.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(ProbioAbunTool, self).run()
        self.run_probio_anno()
        self.end()

    def run_probio_anno(self):
        #out_file = self.output_dir + "/probio_abun.xls"
        probio_file = self.option("probio_anno").prop["path"]
        gene_profile = self.option("reads_profile_table").prop["path"]
        cmd = self.soft + ' {} -i {} -p {} -o {}'.format(self.script, probio_file, gene_profile,self.output_dir)
        self.logger.info(cmd)
        command = self.add_command('probio_abun', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('run probio_abun succeed')
        else:
            self.set_error("probio_abun failed", code="31204401")
            self.set_error('probio_abun failed', code="31204402")

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "probio_abun",
            "type": "tool",
            "name": "annotation.probio_abun",
            "instant": True,
            "options": dict(
                probio_anno="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/metag_v2/nr/gene_probio_anno.xls",
                reads_profile_table="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/annotation/gene_profile.reads_number.total.txt"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
