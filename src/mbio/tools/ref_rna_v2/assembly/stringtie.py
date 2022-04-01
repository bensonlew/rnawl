# -*- coding: utf-8 -*-
# __author__ = "qinjincheng"

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool


class StringtieAgent(Agent):
    def __init__(self, parent):
        super(StringtieAgent, self).__init__(parent)
        options = [
            {"name": "sample_bam", "type": "infile", "format": "align.bwa.bam"},
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "cpu", "type": "int", "default": 8},
            {"name": "fr_stranded", "type": "string", "default": "fr-unstranded"},
            {"name": "strand_direct", "type": "string", "default": "none"},
            {"name": "sample_gtf", "type": "outfile", "format": "gene_structure.gtf"},
            {"name": "min_coverage", "type": "int", "default": 3},
            {"name": "min_read", "type": "int", "default": 5},
            {"name": "verbose", "type": "bool", "default": True},
        ]
        self.add_option(options)
        self._memory_increase_step = 20

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format(int(os.path.getsize(self.option("ref_gtf").path) / 1024.0 ** 3 * 8 + 16))

    def end(self):
        super(StringtieAgent, self).end()


class StringtieTool(Tool):
    def __init__(self, config):
        super(StringtieTool, self).__init__(config)
        self.program = {
            "stringtie": "bioinfo/ref_rna_v2/miniconda2/bin/stringtie"
        }
        sample_name = os.path.basename(self.option("sample_bam").path)[:-4]
        self.file = {
            "gtf": os.path.join(self.output_dir, "{}.gtf".format(sample_name)),
            "gene_abund": os.path.join(self.work_dir, "{}.gene_abund.out".format(sample_name))
        }

    def run(self):
        super(StringtieTool, self).run()
        self.run_stringtie()
        self.end()

    def run_stringtie(self):
        cmd = "{} {}".format(self.program["stringtie"], self.option("sample_bam").path)
        cmd += " -G {}".format(self.option("ref_gtf").path)
        if self.option("fr_stranded") != "fr-unstranded":
            if self.option("strand_direct") == "firststrand":
                cmd += " --rf"
            elif self.option("strand_direct") == "secondstrand":
                cmd += " --fr"
        cmd += " -o {}".format(self.file["gtf"])
        cmd += " -j {}".format(self.option("min_coverage"))
        cmd += " -c {}".format(self.option("min_read"))
        if self.option("verbose"):
            cmd += " -v"
        cmd += " -p {}".format(self.option("cpu"))
        cmd += " -A {}".format(self.file["gene_abund"])
        cmd += " -b {}".format(self.work_dir)
        command = self.add_command("run_stringtie", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.set_output()
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running stringtie")
                self.set_output()
        else:
            self.set_error("fail to run stringtie")

    def set_output(self):
        self.option("sample_gtf").set_path(self.file["gtf"])


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "stringtie_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "ref_rna_v2.assembly.stringtie",
            "instant": False,
            "options": {
                "sample_bam": "/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/RnaseqMapping/output"
                              "/bam/TR1_3.bam",
                "ref_gtf": "/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/FileCheck"
                           "/Arabidopsis_thaliana.TAIR10.43.gtf",
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
