# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
__author__ = 'fengyitong'


class SrnaAnnotationAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(SrnaAnnotationAgent, self).__init__(parent)
        options = [
            {'name': 'predict_fa', 'type': 'string', 'default': ''},
            {"name": "rfam", "type": 'string', 'default': ''},
            {'name': 'evalue', 'type': 'string', 'default': '0.00001'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 5
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
        super(SrnaAnnotationAgent, self).end()

class SrnaAnnotationTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(SrnaAnnotationTool, self).__init__(config)
        self.python_path = '/program/Python/bin/python'
        if self.option('rfam'):
            self.srna_annot = self.config.PACKAGE_DIR + "/prok_rna/sRNA_annotation_v3.py"
        else:
            self.srna_annot = self.config.PACKAGE_DIR + "/prok_rna/sRNA_annotation.py"

    def srnaAnnotation(self):
        cmd = '{} {} '.format(self.python_path, self.srna_annot)
        cmd += '{} {} '.format(self.option("predict_fa"), self.option("evalue"))
        if self.option('rfam'):
            cmd += '{}'.format(self.option('rfam'))
        cmd_name = 'srna_annotation'
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
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004101")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004102")

    def set_output(self):
        all_files = os.listdir(self.work_dir)
        all_files = [self.work_dir + '/' + each for each in all_files ]
        for each in all_files:
            if each.endswith('.xls') or each.endswith('_merge') or each.endswith('list'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(SrnaAnnotationTool, self).run()
        self.srnaAnnotation()
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
            "id": "SrnaAnnotation" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "prok_rna.sRNA_annotation",
            "instant": False,
            "options": dict(
                predict_fa = test_dir + "/" + "genome.predicted_RNA.fa",
                evalue = "0.001",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()