# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
__author__ = 'fengyitong'


class SrnaFoldAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(SrnaFoldAgent, self).__init__(parent)
        options = [
            {'name': 'predict_fa', 'type': 'string', 'default': ''},
        ]
        self.add_option(options)

    def check_options(self):
        pass

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
        super(SrnaFoldAgent, self).end()

class SrnaFoldTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(SrnaFoldTool, self).__init__(config)
        self.python_path = '/program/Python/bin/python'
        self.srna_annot = self.config.PACKAGE_DIR + "/prok_rna/sRNA_fold.py"

    def SrnaFold(self):
        cmd = '{} {} {}'.format(self.python_path, self.srna_annot, self.option("predict_fa"))
        cmd_name = 'srna_fold'
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
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004201")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004202")

    def set_output(self):
        all_files = os.listdir(self.work_dir + '/' + 'RNAfold_pdf')
        all_files = [self.work_dir + '/RNAfold_pdf/' + each for each in all_files ]
        for each in all_files:
            if each.endswith('.pdf') or each.endswith('txt') or each.endswith('.ps'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
        # os.link(os.path.join(self.work_dir, 'RNAfold.str.txt'), os.path.join(self.output_dir, 'RNAfold.str.txt'))
        os.link(os.path.join(self.work_dir, 'RNAfold.str.txt'), os.path.join(self.output_dir, 'RNAfold.str'))

    def run(self):
        super(SrnaFoldTool, self).run()
        self.SrnaFold()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/sRNA_annotation'
        data = {
            "id": "SrnaFold_" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "prok_rna.sRNA_fold",
            "instant": False,
            "options": dict(
                predict_fa = test_dir + "/" + "genome.predicted_RNA.fa",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()