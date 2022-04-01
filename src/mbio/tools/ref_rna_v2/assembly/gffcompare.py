# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool


class GffcompareAgent(Agent):

    def __init__(self, parent):
        super(GffcompareAgent, self).__init__(parent)
        options = [
            {"name": "merged_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "cuff_gtf", "type": "outfile", "format": "gene_structure.gtf"},
            {"name": "tmap", "type": "outfile", "format": "assembly.tmap"},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "4G"

    def end(self):
        super(GffcompareAgent, self).end()


class GffcompareTool(Tool):
    def __init__(self, config):
        super(GffcompareTool, self).__init__(config)
        self.program = {
            "gffcompare": "bioinfo/rna/gffcompare-0.9.8.Linux_x86_64/gffcompare"
        }
        self.file = {
            'input_gtf': os.path.join(self.work_dir, 'input.gtf'),
            'tmap': os.path.join(self.work_dir, 'gffcmp.input.gtf.tmap')
        }
        if not os.path.isfile(self.file['input_gtf']):
            os.link(self.option('merged_gtf').path, self.file['input_gtf'])

    def run(self):
        super(GffcompareTool, self).run()
        self.run_gffcompare()
        self.end()

    def run_gffcompare(self):
        cmd = '{} -r {} {}'.format(self.program['gffcompare'], self.option('ref_gtf').path, self.file['input_gtf'])
        command = self.add_command("run_gffcompare", cmd)
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
                self.logger.info("succeed in running gffcompare")
                self.set_output()
        else:
            self.set_error("fail to run gffcompare")

    def set_output(self):
        link_name = os.path.join(self.output_dir, os.path.basename(self.file['tmap']))
        if os.path.isfile(link_name):
            os.remove(link_name)
        os.link(self.file['tmap'], link_name)
        self.option('tmap').set_path(link_name)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "gffcompare_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "ref_rna_v2.assembly.gffcompare",
            "instant": False,
            "options": {
                "merged_gtf": "/mnt/ilustre/users/sanger-dev/workspace/20200324/Single_stringtie_merge_6944_9677"
                              "/StringtieMerge/output/out.gtf",
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
