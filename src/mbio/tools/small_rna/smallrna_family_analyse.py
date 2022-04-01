# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
__author__ = 'fengyitong'


class SmallrnaFamilyAnalyseAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(SmallrnaFamilyAnalyseAgent, self).__init__(parent)
        options = [
            {'name': 'family', 'type': 'string', 'default': self.config.SOFTWARE_DIR + '/database/mirbase/family_analyse/miRfamily.dat'},
            {'name': 'pre', 'type': 'string', 'default': self.config.SOFTWARE_DIR + '/database/mirbase/family_analyse/miR_pre_mature.dat'},
            {'name': 'mir', 'type': 'string'},
            {'name': 'matfa', 'type': 'string'},
            {'name': 'novofa', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('5')

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
        super(SmallrnaFamilyAnalyseAgent, self).end()

class SmallrnaFamilyAnalyseTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(SmallrnaFamilyAnalyseTool, self).__init__(config)
        self.python_path = 'program/Python/bin/python'
        self.FamilyAnalyse = self.config.PACKAGE_DIR + "/small_rna/smallrna_family_analyse.py"

    def run_FamilyAnalyse(self):
        cmd = '{} {} '.format(self.python_path, self.FamilyAnalyse)
        cmd += '-{} {} '.format("fam", self.option("family"))
        cmd += '-{} {} '.format("pre", self.option("pre"))
        cmd += '-{} {} '.format("mir", self.option("mir"))
        cmd += '-{} {} '.format("matfa", self.option("matfa"))
        cmd += '-{} {} '.format("novofa", self.option("novofa"))
        cmd_name = 'family_analyse'
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
                each = os.path.join(self.work_dir, fname)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(SmallrnaFamilyAnalyseTool, self).run()
        self.run_FamilyAnalyse()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/smallrna/family'
        data = {
            "id": "FamilyAnalyse" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.smallrna_family_analyse",
            "instant": False,
            "options": dict(
                # mir=test_dir + "/" + "miRNAs_expressed.xls",
                # matfa=test_dir + "/" + "mature_mmu.dna.fa",
                # novofa=test_dir + "/" + "novo_mature_seq.fa",
                family="/mnt/ilustre/users/sanger-dev/app/database/PmiREN/family_analyse/miRfamily.dat",
                pre="/mnt/ilustre/users/sanger-dev/app/database/PmiREN/family_analyse/miR_pre_mature.dat",
                mir="/mnt/ilustre/users/sanger-dev/workspace/20200508/Single_KnownMirna_2212/KnownMirna/output/known_mirna_count.xls",
                matfa="/mnt/ilustre/users/sanger-dev/workspace/20200508/Single_KnownMirna_2212/KnownMirna/output/mature.fa",
                novofa="/mnt/ilustre/users/sanger-dev/workspace/20200508/Single_NovelMirna_2983/NovelMirna/output/novel_mature_seq.fa"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()