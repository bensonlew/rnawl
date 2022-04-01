# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from Bio import SeqIO
# import pandas as pd
__author__ = '刘彬旭'


class NormalizeFastqAgent(Agent):
    """
    spades assemble
    """
    def __init__(self, parent):
        super(NormalizeFastqAgent, self).__init__(parent)
        options = [
            {'type': 'string', 'name': 'mode', 'default': 'paired'},
            {'type': 'infile', 'name': 'l', 'format': 'denovo_rna_v2.common'},
            {'type': 'infile', 'name': 'r', 'format': 'denovo_rna_v2.common'},
            {'type': 'infile', 'name': 's', 'format': 'denovo_rna_v2.common'},
            {'type': 'string', 'name': 'o', 'default': ''},
            {'type': 'string', 'name': 'soft', 'default': 'ORNA'},
            {'type': 'int', 'name': 't', 'default': 20},
            {'type': 'string', 'name': 'k', 'default': "auto"},
            {'type': 'string', 'name': 'c', 'default': "auto"},
            {'type': 'int', 'name': 'Q', 'default': 'auto'},
            {'type': 'string', 'name': 'n', 'default': 'NODE'},
        ]
        self.add_option(options)

    def check_options(self):
        self._memory_increase_step = 20
        pass


    def set_resource(self):
        file_size = float(os.path.getsize(self.option('l').prop['path'])) / 1024 / 1024 /1024
        if self.option('soft') == "Trinity":
            mem = "200G"
        else:
            mem = "50G"
        self._memory = "{}G".format(mem)
        self._cpu = self.option("t") + 1

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
        super(NormalizeFastqAgent, self).end()

class NormalizeFastqTool(Tool):
    """
    spades assemble
    """
    def __init__(self, config):
        super(NormalizeFastqTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.orna = 'bioinfo/denovo_rna_v2/miniconda2/pkgs/orna-2.0-he52c88d_0/bin/ORNA'
        self.trinity = 'bioinfo/denovo_rna_v2/trinityrnaseq-2.8.5/util/insilico_read_normalization.pl'
        self.jellyfish = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/jellyfish-2.3.0/bin/'
        self._LD_LIBRARY_PATH = software_dir + "/bioinfo/denovo_rna_v2/miniconda2/lib"
        self.set_environ(LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.set_environ(PATH=self.jellyfish)

    def run_normalization_orna(self):        
        if self.option("mode") == "paired":
            cmd = "{} -pair1  {} -pair2 {} -output {} -nb-cores {} -type fastq".format(
                self.orna, 
                self.option("l").prop['path'], 
                self.option("r").prop['path'], 
                self.output_dir + "/normalize",
                self.option("t")
            )
        else:
            cmd = "{} -input {}  -output {} -nb-cores {} -type fastq".format(
                self.orna, 
                self.option("l").prop['path'], 
                self.output_dir + "/normalize",
                self.option("t")
            ) 

        cmd_name = 'orna'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        else:
            self.set_error("{} Failed and returned None, we will try it again?".format(cmd_name))

    def run_normalization_trinity(self):
        if self.option("mode") == "paired":
            cmd = "{} --seqType fq --JM 190G  --max_cov 200 --min_cov 5 --max_CV 10000  --left {} --right {} --output {} --CPU {} --pairs_together --PARALLEL_STATS".format(
                self.trinity, 
                self.option("l").prop['path'], 
                self.option("r").prop['path'], 
                self.work_dir + "/normalize",
                self.option("t")
            )
        else:
            cmd = "{} --seqType fq --JM 190G  --max_cov 200 --min_cov 5 --max_CV 10000  --single {} --output {} --CPU {} --PARALLEL_STATS".format(
                self.trinity, 
                self.option("l").prop['path'], 
                self.work_dir + "/normalize",
                self.option("t")
            )

        cmd_name = 'trinity_norm'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        else:
            self.set_error("{} Failed and returned None, we will try it again?".format(cmd_name))
    def set_output(self):
        if self.option("soft") == "Trinity":
            os.link(
                os.path.join(self.work_dir, "normalize", "left.norm.fq"),
                os.path.join(self.output_dir, "normalize_1.fq")
            )
            os.link(
                os.path.join(self.work_dir, "normalize", "right.norm.fq"),
                os.path.join(self.output_dir, "normalize_2.fq")
            )
        pass

    def run(self):
        super(NormalizeFastqTool, self).run()
        if self.option("soft") == "ORNA":
            self.run_normalization_orna()
        elif self.option("soft") == "Trinity":
            self.run_normalization_trinity()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/ref_denovo/'
        data = {
            "id": "spades" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.spades",
            "instant": False,
            "options": dict(
                l=test_dir + "C_0_5R_1_S211_L004_R1_001.fastq",
                r=test_dir + "C_0_5R_1_S211_L004_R2_001.fastq",
                n='sample1'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
