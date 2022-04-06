# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import json
import os
import shutil
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class ProteinPictureAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(ProteinPictureAgent, self).__init__(parent)
        options = [
            {'name': 'picture', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'if_black', 'type': 'string', 'default': 'no'},
            {'name': 'protein_info_table', 'type': 'outfile', 'format': 'medical_transcriptome.common'},

        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'
    def end(self):
        super(ProteinPictureAgent, self).end()


class ProteinPictureTool(Tool):
    def __init__(self, config):
        super(ProteinPictureTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'picture': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/get_protein_infomation_from_picture.py'),
        }
        self.file = {
            'protein_info_table': os.path.join(self.output_dir, 'Protein_information.xls'),
        }

    def run(self):
        super(ProteinPictureTool, self).run()
        self.protein_info()
        self.set_output()
        self.end()

    def protein_info(self):
        cmd = '{} {}'.format(self.program['python'], self.script['picture'])
        cmd += ' -picture_path {}'.format(self.option('picture').path)
        cmd += ' -if_black {}'.format(self.option('if_black'))
        cmd += ' -output {}'.format(self.output_dir)
        cmd_name = 'run_protein_info'
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
        self.option('protein_info_table').set_path(self.file['protein_info_table'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'protein_picture_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.protein_picture',
            'instant': False,
            'options': {
                'picture': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/protein_picture/Dl09YYJ_bjhb1_1580_Result.PNG',
                'if_black': 'no',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
