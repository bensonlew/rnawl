# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
__author__ = 'fengyitong'


class SrnaTargetAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(SrnaTargetAgent, self).__init__(parent)
        options = [
            {'name': 'predict_fa', 'type': 'string', 'default': ''},
            {'name': 'genome_bed', 'type': 'string', 'default': ''},
            {'name': 'genome_fa', 'type': 'string', 'default': ''},
        ]
        self.add_option(options)
        self._memory_increase_step = 50

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = "{}G".format('80')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".xls", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(SrnaTargetAgent, self).end()

class SrnaTargetTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(SrnaTargetTool, self).__init__(config)
        self.python_path = '/miniconda2/bin/python'
        self.srna_annot = self.config.PACKAGE_DIR + "/prok_rna/sRNA_target.py"

    def SrnaTarget(self):
        cmd = '{} {} '.format(self.python_path, self.srna_annot)
        cmd += '{} {} {}'.format(self.option("predict_fa"), self.option("genome_bed"), self.option("genome_fa"))
        cmd_name = 'srna_target'
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
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004301")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004302")

    def set_output(self):
        all_files = os.listdir(self.work_dir)
        all_files = [self.work_dir + '/' + each for each in all_files ]
        for each in all_files:
            if each.endswith('merge') or os.path.basename(each).startswith('combine'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(SrnaTargetTool, self).run()
        self.SrnaTarget()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/target_test'
        data = {
            "id": "SrnaTarget_" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "prok_rna.sRNA_target",
            "instant": False,
            "options": dict(
                predict_fa = test_dir + "/" + "genome.predicted_RNA.fa",
                genome_bed = test_dir + "/" + "genome.gene.bed",
                genome_fa = test_dir + "/" + "GCF_000009345.1_ASM934v1_genomic.fna",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()