# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import pandas as pd
import numpy as np
from sklearn import decomposition, preprocessing


class PlantNameTranslationAgent(Agent):
    """
    Used for convert RNA-seq expression units convert, count/CPM/RPM/TPM/FPKM/RPKM are implemented.

    RPM/CPM: Reads/Counts of exon model per million mapped reads
    RPM/CPM=Total exon reads/ Mapped reads(Millions)

    RPKM/FPKM: Reads/Fragments Per Kilobase of exon model per Million mapped reads
    RPKM/FPKM=Total exon reads/[Mapped reads(Millions)*Exon length(Kb)]

    TPM is like RPKM and FPKM, except the order of operation is switched.
    """
    def __init__(self, parent):
        super(PlantNameTranslationAgent, self).__init__(parent)
        options = [
            {'name': 'name_type', 'type': 'string'},
            {'name': 'plant_list', 'type': 'string'},
            {'name': 'translation_file', 'type': 'outfile', 'format': 'medical_transcriptome.common'}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "3G"

    def end(self):
        super(PlantNameTranslationAgent, self).end()


class PlantNameTranslationTool(Tool):
    def __init__(self, config):
        super(PlantNameTranslationTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'translation': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/plant_name_translation.py'),
        }
        self.file = {
            'output_file': os.path.join(self.output_dir, 'translation.txt')
        }
    def run(self):
        """
        运行
        :return:
        """
        super(PlantNameTranslationTool, self).run()
        self.run_translation()
        self.set_output()
        self.end()

    def run_translation(self):
        cmd = '{} {}'.format(self.program['python'], self.script['translation'])
        cmd += ' -n {}'.format(self.option('name_type'))
        cmd += ' -p "{}"'.format(str(self.option('plant_list')))
        cmd += ' -o {}'.format(self.file['output_file'])
        cmd_name = 'run_plant'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704403")

    def set_output(self):
        self.option('translation_file').set_path(self.file['output_file'])

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "heatmap_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.heatmap",
            "options": dict(
                otutable="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/cluster/otu_table.xls",
                log_change="None",
                scale='zscore'

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
