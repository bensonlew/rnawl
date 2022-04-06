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


class StrCasecontrolAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(StrCasecontrolAgent, self).__init__(parent)
        options = [
            {'name': 'manifest', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'merge_json', 'type': 'infile', 'format': 'ref_rna_v2.common'}
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
        super(StrCasecontrolAgent, self).end()


class StrCasecontrolTool(Tool):
    def __init__(self, config):
        super(StrCasecontrolTool, self).__init__(config)
        self.program = {
            'ExpansionHunterDenovo': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/script/ExpansionHunterDenovo/ExpansionHunterDenovo-v0.9.0-linux_x86_64/bin/ExpansionHunterDenovo',
            'casecontrol': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/script/ExpansionHunterDenovo/ExpansionHunterDenovo-v0.9.0-linux_x86_64/scripts/casecontrol.py',
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'combat': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/batch/combat.r'),
            'removebtacheffect': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/batch/removebatcheffect.r')
        }
        self.file = {
            'output_locus': os.path.join(self.output_dir, 'str.casecontrol_locus.tsv'),
            'output_motif': os.path.join(self.output_dir, 'str.casecontrol_motif.tsv')
        }

    def run(self):
        super(StrCasecontrolTool, self).run()
        self.run_str_casecontrol_locus()
        self.run_str_casecontrol_motif()
        # self.set_output()
        self.end()

    def run_str_casecontrol_locus(self):
        cmd = '{} {} '.format(self.program['python'], self.program['casecontrol'])
        cmd += '{} '.format('locus')
        cmd += '--manifest {} '.format(self.option('manifest').path)
        cmd += '--multisample-profile {} '.format(self.option('merge_json').path)
        cmd += '--output {} '.format(self.file['output_locus'])
        cmd_name = 'run_casecontrol_locus'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    def run_str_casecontrol_motif(self):
        cmd = '{} {} '.format(self.program['python'], self.program['casecontrol'])
        cmd += '{} '.format('motif')
        cmd += '--manifest {} '.format(self.option('manifest').path)
        cmd += '--multisample-profile {} '.format(self.option('merge_json').path)
        cmd += '--output {} '.format(self.file['output_motif'])
        cmd_name = 'run_casecontrol_motif'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")


    def set_output(self):
        self.option('merge_json').set_path(self.file['merge_json'])



class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'STR{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.str_casecontrol',
            'instant': False,
            'options': {
                'manifest': '/mnt/ilustre/users/sanger-dev/workspace/20210121/Single_Str_predict_7590_3534/Str/manifest.tsv',
                'merge_json': '/mnt/ilustre/users/sanger-dev/workspace/20210121/Single_Str_predict_7590_3534/Str/StrMerge/output/str.multisample_profile.json',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
