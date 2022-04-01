# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import ConfigParser
import unittest
__author__ = 'fengyitong'


class MirnaTarget2Agent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(MirnaTarget2Agent, self).__init__(parent)
        options = [
            {'name': 'query', 'type': 'string'},
            {'name': 'target', 'type': 'string'},
            {'name': 'method', 'type': 'string'},
            {'name': 'spiece', 'type': 'string'},
            {'name': 'target_species', 'type': 'string', 'default': 'Homo_sapiens'},
            {"name": "miranda_score", "type": "string", "default": "160"},
            {"name": "miranda_energy", "type": "string", "default": "-20"},
            {"name": "miranda_strict", "type": "string", "default": "on"},
            {"name": "rnahybird_num", "type": "string", "default": "100"},
            {"name": "rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "rnahybrid_num", "type": "string", "default": "100"},
            {"name": "rnahybrid_energy", "type": "string", "default": "-20"},
            {"name": "rnahybrid_pvalue", "type": "string", "default": "0.01"},
            {"name": "ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "targetfinder_score", "type": "string", "default": "4"},

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 6
        self._memory = "{}G".format('60')

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
        super(MirnaTarget2Agent, self).end()

class MirnaTarget2Tool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(MirnaTarget2Tool, self).__init__(config)
        self.python_path = 'program/Python/bin/python'
        self.MirnaTarget = self.config.PACKAGE_DIR + "/small_rna/mirna_target2.py"
        self.perl = '/program/perl-5.24.0/bin/'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        self.ssearch = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/fasta-35.4.12/bin'
        self.set_environ(PATH=self.ssearch)
        self.software_dict = dict(
        miranda = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/miRanda-3.3a/bin/miranda",
        pita = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/pita_v6/pita_prediction.pl",
        rnahybrid = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/RNAhybrid-2.1.2/src/RNAhybrid",
        psrobot = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/psRobot_v1.2/psRobot_tar",
        targetfinder = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/targetfinder/targetfinder_threads.pl",
        targetscan = self.config.SOFTWARE_DIR + "/database/targetscan/Predicted_Targets_Info.default_predictions.txt"
        )
        self.target_context = self.config.SOFTWARE_DIR + "/database/targetscan/Predicted_Targets_Context_Scores.default_predictions.txt"

    def run_MirnaTarget(self):
        cmd = '{} {} '.format(self.python_path, self.MirnaTarget)
        cmd += '-q {} '.format(self.option('query'))
        cmd += '-t {} '.format(self.option('target'))
        cmd += '-m {} '.format(self.option('method').lower())
        cmd += '-miranda_s {} '.format(self.option('miranda_score'))
        cmd += '-miranda_e {} '.format(self.option('miranda_energy'))
        cmd += '-miranda_strict {} '.format(self.option('miranda_strict'))
        cmd += '-rnahybrid_b {} '.format(self.option('rnahybrid_num'))
        cmd += '-rnahybrid_e {} '.format(self.option('rnahybrid_energy'))
        cmd += '-rnahybrid_p {} '.format(self.option('rnahybrid_pvalue'))
        cmd += '-psrobot_ts {} '.format(self.option('ps_robot_score'))
        cmd += '-targetfinder_c {} '.format(self.option('targetfinder_score'))
        cmd += '-spiece {} '.format(self.option('spiece'))
        cmd += '-{} {} '.format(self.option('method').lower(),
                                self.software_dict[self.option('method').lower()])

        if self.option('method') == 'targetscan':
            spe_tax = {
                "Mus_musculus": 10090,
                "Rattus_norvegicus": 10116,
                "Monodelphis_domestica": 13616,
                "Xenopus_tropicalis": 8364,
                "Gallus_gallus": 9031,
                "Macaca_mulatta": 9544,
                "Pan_troglodytes": 9598,
                "Homo_sapiens":	9606,
                "Canis_lupus_familiaris": 9615,
                "Bos_taurus": 9913
            }

            if not self.option('target_species') in spe_tax:
                raise Exception("species un supported must in {}".format(spe_tax.keys()))
            cmd += '-{} {} '.format("target_species", spe_tax[self.option('target_species')])
            cmd += '-{} {} '.format("target_context", self.target_context)
        cmd_name = 'mirnatarget'
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
            if each.endswith('merge_out') or each.endswith("detail.txt.gz"):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(MirnaTarget2Tool, self).run()
        self.run_MirnaTarget()
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
            "id": "MirnaTarget2_psrobot" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.mirna_target2",
            "instant": False,
            "options": dict(
                # target='/mnt/ilustre/users/sanger-dev/sg-users/machao/miRNA_target_test/animal_data/Rnor_6.0.91.3utr.fa',
                target='/mnt/ilustre/users/sanger-dev/sg-users/machao/miRNA_target_test/plant_data/ref_genome.gtf.exon.fa',
                # query='/mnt/ilustre/users/sanger-dev/sg-users/machao/miRNA_target_test/animal_data/all_known_DEM.fa',
                query='/mnt/ilustre/users/sanger-dev/sg-users/machao/miRNA_target_test/plant_data/mature_ptc.dna.fa',
                spiece='plant',
                # spiece='plant',
                method='psrobot'
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
