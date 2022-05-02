# -*- coding: utf-8 -*-
# __author__ = "qinjincheng"

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool


class StringtieMergeAgent(Agent):
    def __init__(self, parent):
        super(StringtieMergeAgent, self).__init__(parent)
        options = [
            {"name": "gtf_list_fp", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "ref_fa", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "cpu", "type": "int", "default": 10},
            {"name": "merged_gtf", "type": "outfile", "format": "gene_structure.gtf"},
            {"name": "min_iso", "type": "float", "default": 0.1},
            {"name": "min_cov", "type": "int", "default": 5},
            {"name": "min_tpm", "type": "int", "default": 1},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "8G"

    def end(self):
        super(StringtieMergeAgent, self).end()


class StringtieMergeTool(Tool):
    def __init__(self, config):
        super(StringtieMergeTool, self).__init__(config)
        self.program = {
            "stringtie": "miniconda2/bin/stringtie"
        }
        self.file = {
            "out_gtf": os.path.join(self.output_dir, "out.gtf"),
        }

    def run(self):
        super(StringtieMergeTool, self).run()
        self.run_stringtie_merge()
        self.end()

    def run_stringtie_merge(self):
        cmd = "{} --merge".format(self.program["stringtie"])
        cmd += " -G {}".format(self.option("ref_gtf").path)
        cmd += " -o {}".format(self.file["out_gtf"])
        cmd += " -c {}".format(self.option("min_cov"))
        cmd += " -T {}".format(self.option("min_tpm"))
        cmd += " -f {}".format(self.option("min_iso"))
        cmd += ' {}'.format(self.option('gtf_list_fp').path)
        command = self.add_command("run_stringtie_merge", cmd)
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
        self.option("merged_gtf").set_path(self.file["out_gtf"])


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "stringtie_merge_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "ref_rna_v2.assembly.stringtie_merge",
            "instant": False,
            "options": {
                "gtf_list_fp": "/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/RefrnaAssemble"
                               "/assembly_gtf.txt",
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
