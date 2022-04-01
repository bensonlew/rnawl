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


class StrExpansionhunterAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(StrExpansionhunterAgent, self).__init__(parent)
        options = [
            {'name': 'bam', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'bam_name', 'type': 'string'},
            {'name': 'variant-catalog', 'type': 'infile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)


    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 5
        self._memory = '50G'

    def end(self):
        super(StrExpansionhunterAgent, self).end()


class StrExpansionhunterTool(Tool):
    def __init__(self, config):
        super(StrExpansionhunterTool, self).__init__(config)
        # software_dir = self.config.SOFTWARE_DIR
        # self.gcc = software_dir + '/gcc/5.1.0/bin'
        # self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        # self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'samtools': '/bioinfo/align/samtools-1.3.1/samtools',
            'ExpansionHunter': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/script/ExpansionHunter/ExpansionHunter-v4.0.2-linux_x86_64/bin/ExpansionHunter'
        }
        self.script = {
            'combat': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/batch/combat.r'),
            'removebtacheffect': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/batch/removebatcheffect.r')
        }

        self.file = {
            'output_prefix': os.path.join(self.output_dir, self.option('bam_name')),
        }

    def run(self):
        super(StrExpansionhunterTool, self).run()
        if not os.path.isfile('{}.{}'.format(self.option('bam').path, 'bai')):
            self.run_samtools_index()
        self.run_str()
        self.set_output()
        self.end()

    def run_samtools_index(self):
        cmd = '{} {} -b '.format(self.program['samtools'], 'index')
        cmd += '{}'.format(self.option('bam').path)
        cmd_name = 'run_samtools_index'
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
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704307")

    def run_str(self):
        cmd = '{} '.format(self.program['ExpansionHunter'])
        cmd += '--reads {} '.format(self.option('bam').path)
        cmd += '--reference {} '.format(self.option('ref').path)
        cmd += '--output-prefix {} '.format(self.file['output_prefix'])
        cmd += '--variant-catalog {} '.format(self.option('variant-catalog').path)
        cmd_name = 'run_expansionhunter'
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
        # self.option('str_json').set_path(self.file['output_json'])
        pass



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
            'name': 'tool_lab.str_expansionhunter',
            'instant': False,
            'options': {
                'bam': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/data/sample2.bam',
                'ref': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/SL4.0_ITAG4.0/dna/S_lycopersicum_chromosomes.4.00.fa',
                'variant-catalog': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/data/variant_catalog/variant_catalog_ssr_exon.json'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
