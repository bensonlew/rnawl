# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class CmscanAgent(Agent):
    '''
    last_modify: 2019.04.19
    '''

    def __init__(self, parent):
        super(CmscanAgent, self).__init__(parent)
        options = [
            {'name': 'seqfile', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'fmt', 'type': 'int', 'default': 2},
            {'name': 'cpu', 'type': 'int', 'default': 8},
            {'name': 'tsv', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('cmscan')
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.cmscan.start()
        self.step.update()

    def step_finish(self):
        self.step.cmscan.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('seqfile').is_set:
            self.logger.debug('{} = {}'.format('seqfile', self.option('seqfile').prop['path']))
            self.infile_size = os.path.getsize(self.option('seqfile').prop['path'])
        else:
            raise OptionError('query sequences in FASTA format must be provided')
        self.logger.debug('{} = {}'.format('fmt', self.option('fmt')))
        self.logger.debug('{} = {}'.format('cpu', self.option('cpu')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('cpu')
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 8 + 5))

    def end(self):
        super(CmscanAgent, self).end()


class CmscanTool(Tool):
    def __init__(self, config):
        super(CmscanTool, self).__init__(config)
        self.cmscan = 'bioinfo/lnc_rna/miniconda2/bin/cmscan'
        self.python = 'program/Python/bin/python'
        self.scantbl2tsv = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/scantbl2tsv.py')
        self.tblout = os.path.join(self.work_dir, 'cmscan.tblout')
        self.Rfam_clanin = os.path.join(self.config.SOFTWARE_DIR, 'database/lnc_rna/Rfam/Rfam.clanin')
        self.Rfam_with_desc = os.path.join(self.config.SOFTWARE_DIR, 'database/lnc_rna/Rfam/Rfam_with_desc.cm')
        self.tsv = os.path.join(self.output_dir, 'cmscan.tsv')

    def run(self):
        super(CmscanTool, self).run()
        self.run_cmscan()
        self.run_scantbl2tsv()
        self.set_output()
        self.end()

    def run_cmscan(self):
        cmd = '{} --rfam --cut_ga --nohmmonly --notextw --oskip'.format(self.cmscan)
        cmd += ' --fmt {}'.format(self.option('fmt'))
        cmd += ' -o {} '.format(os.path.join(self.work_dir, 'cmscan.stdout'))
        cmd += ' --tblout {} '.format(self.tblout)
        cmd += ' --clanin {} '.format(self.Rfam_clanin)
        cmd += ' --cpu {} '.format(self.option('cpu'))
        cmd += ' {} {}'.format(self.Rfam_with_desc, self.option('seqfile').prop['path'])
        cmd_name = 'run_cmscan'
        self.run_code(cmd_name, cmd)

    def run_scantbl2tsv(self):
        cmd = '{} {}'.format(self.python, self.scantbl2tsv)
        cmd += ' -i {}'.format(self.tblout)
        cmd += ' -o {}'.format(self.tsv)
        cmd_name = 'run_scantbl2tsv'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run {}, abord'.format(n))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('tsv').set_path(self.tsv)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'cmscan_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.cmscan',
            'instant': False,
            'options': {
                'seqfile': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/FilterByExpress/output/filtered_file/all_lncrna.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
