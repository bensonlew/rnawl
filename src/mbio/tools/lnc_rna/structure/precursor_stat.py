# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class PrecursorStatAgent(Agent):
    '''
    last_modify: 2019.03.27
    '''
    def __init__(self, parent):
        super(PrecursorStatAgent, self).__init__(parent)
        options = [
            {'name': 'tblout', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'species', 'type': 'string', 'default': ''},
            {'name': 'known_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'novel_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 't2g', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'tabular', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
        ]
        self.add_option(options)
        self.step.add_steps('precursor_stat')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.precursor_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.precursor_stat.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        super(PrecursorStatAgent, self).end()

class PrecursorStatTool(Tool):
    def __init__(self, config):
        super(PrecursorStatTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.hairpin_fa = os.path.join(
            self.config.SOFTWARE_DIR,
            'database/mirbase/species/{}/{}.hairpin.fa'.format(self.option('species'), self.option('species'))
        )
        self.precursor_statistics_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/precursor_statistics.py')
        self.tabular = os.path.join(self.work_dir, 'miRNA_precursor.tabular')

    def run(self):
        super(PrecursorStatTool, self).run()
        self.run_precursor_stat()
        self.set_output()
        self.end()

    def run_precursor_stat(self):
        cmd = '{} {}'.format(self.python, self.precursor_statistics_py)
        cmd += ' -i {}'.format(self.option('tblout').path)
        cmd += ' -k {}'.format(self.option('known_list').path)
        cmd += ' -n {}'.format(self.option('novel_list').path)
        cmd += ' -t {}'.format(self.option('t2g').path)
        cmd += ' -p {}'.format(self.hairpin_fa)
        cmd += ' -o {}'.format(self.tabular)
        cmd_name = 'run_precursor_stat'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.tabular
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('tabular').set_path(link_name)
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
            'id': 'precursor_stat_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.precursor_stat',
            'instant': False,
            'options': {
                'tblout': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/precursor_stat/blastn.tblout',
                'species': 'hsa',
                'known_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/precursor_stat/known_lncrna_ids.list',
                'novel_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/precursor_stat/novel_lncrna_ids.list',
                't2g': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/precursor_stat/trans2gene',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()