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


class PiceditAgent(Agent):
    def __init__(self, parent):
        super(PiceditAgent, self).__init__(parent)
        options = [
            {'name': 'picture_dir', 'type': 'string'},
            {'name': 'edit_file', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "30G"

    def end(self):
        super(PiceditAgent, self).end()


class PiceditTool(Tool):
    def __init__(self, config):
        super(PiceditTool, self).__init__(config)
        self.program = {
            'python': 'bioinfo/ref_rna_v3/hisat2/miniconda3/bin/python',
        }
        self.script = {
            'picedit': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/AutoMarker/PicEdit.py'),
        }
        self.module_file = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/AutoMarker/inference_graph/frozen_inference_graph.pb')
        self.font_path = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/AutoMarker/Founder.ttf')
    def run(self):
        """
        运行
        :return:
        """
        super(PiceditTool, self).run()
        self.run_picedit()
        self.set_output()
        self.end()

    def run_picedit(self):
        cmd = '{} {}'.format(self.program['python'], self.script['picedit'])
        cmd += ' -i {}'.format(self.option('picture_dir'))
        cmd += ' -t {}'.format(self.option('edit_file'))
        cmd += ' -p {}'.format(self.module_file)
        cmd += ' -f {}'.format(self.font_path)
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_picedit'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code in [0, -11]:
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
        # self.option('translation_file').set_path(self.file['output_file'])
        pass

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
