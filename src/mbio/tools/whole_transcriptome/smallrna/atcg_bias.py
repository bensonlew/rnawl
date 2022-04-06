# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
__author__ = 'fengyitong'


class AtcgBiasAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(AtcgBiasAgent, self).__init__(parent)
        options = [
            {'name': 'known', 'type': 'string', 'default': '/mnt/ilustre/users/sanger-dev/sg-users/zhuwengen/whole_transcriptome/smallrna/smallrna_family_analyse/known_mature.fa'},
            {'name': 'novel', 'type': 'string', 'default': '/mnt/ilustre/users/sanger-dev/sg-users/zhuwengen/whole_transcriptome/smallrna/smallrna_family_analyse/novel_mature.fa'},
            # {'name': 'min', 'type': 'string', 'default': '18'},
            # {'name': 'max', 'type': 'string', 'default': '32'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('15')

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
        super(AtcgBiasAgent, self).end()

class AtcgBiasTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(AtcgBiasTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.AtcgBias = self.config.PACKAGE_DIR + "/small_rna/atcg_bias.py"

    def run_AtcgBias(self):
        cmd = '{} {} '.format(self.python_path, self.AtcgBias)
        cmd += '-{} {} '.format("known", self.option("known"))
        cmd += '-{} {} '.format("novel", self.option("novel"))
        cmd_name = 'atcg_bias'
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
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def set_output(self):
        all_files = os.listdir(self.work_dir)
        for each in all_files:
            if each.endswith('.xls'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(AtcgBiasTool, self).run()
        self.run_AtcgBias()
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
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/smallrna/bias/'
        data = {
            "id": "whole_transcriptome_smallrna_AtcgBias" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "whole_transcriptome.smallrna.atcg_bias",
            "instant": False,
            "options": dict(
                known="/mnt/ilustre/users/sanger-dev/sg-users/zhuwengen/whole_transcriptome/smallrna/smallrna_family_analyse/known_mature.fa",
                novel="/mnt/ilustre/users/sanger-dev/sg-users/zhuwengen/whole_transcriptome/smallrna/smallrna_family_analyse/novel_mature.fa",
                # min="18",
                # max="32",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()