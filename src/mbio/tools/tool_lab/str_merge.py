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


class StrMergeAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(StrMergeAgent, self).__init__(parent)
        options = [
            {'name': 'manifest', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'merge_json', 'type': 'outfile', 'format': 'ref_rna_v2.common'}
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
        super(StrMergeAgent, self).end()


class StrMergeTool(Tool):
    def __init__(self, config):
        super(StrMergeTool, self).__init__(config)
        self.program = {
            'ExpansionHunterDenovo': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/script/ExpansionHunterDenovo/ExpansionHunterDenovo-v0.9.0-linux_x86_64/bin/ExpansionHunterDenovo'
        }
        self.script = {
            'combat': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/batch/combat.r'),
            'removebtacheffect': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/batch/removebatcheffect.r')
        }
        self.file = {
            'output': os.path.join(self.output_dir, 'str'),
            'merge_json': os.path.join(self.output_dir, 'str.multisample_profile.json')
        }

    def run(self):
        super(StrMergeTool, self).run()
        self.run_str_merge()
        # self.set_output()
        self.end()

    def run_str_merge(self):
        cmd = '{} {} '.format(self.program['ExpansionHunterDenovo'], 'merge')
        cmd += '--manifest {} '.format(self.option('manifest').path)
        cmd += '--reference {} '.format(self.option('ref').path)
        cmd += '--output-prefix {} '.format(self.file['output'])
        cmd_name = 'run_expansionhunter_merge'
        command = self.add_command(cmd_name, cmd, shell=True)
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
            'name': 'tool_lab.str',
            'instant': False,
            'options': {
                'bam': '/mnt/ilustre/users/sanger-dev/workspace/20210113/Single_RnaseqMapping_3917/RnaseqMapping/output/bam/sample2.bam',
                'ref': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/SL4.0_ITAG4.0/dna/S_lycopersicum_chromosomes.4.00.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
