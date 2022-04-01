# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class AnnotFileAgent(Agent):
    '''
    last_modify: 2019.02.22
    '''
    def __init__(self, parent):
        super(AnnotFileAgent, self).__init__(parent)
        options = [
            {'name': 'gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'g2t2p', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'gene2trans', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'i2u', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('annot_file')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.annot_file.start()
        self.step.update()

    def step_end(self):
        self.step.annot_file.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('gtf').is_set:
            self.logger.debug('{} = {}'.format('gtf', self.option('gtf').prop['path']))
        if self.option('g2t2p').is_set:
            self.logger.debug('{} = {}'.format('g2t2p', self.option('g2t2p').prop['path']))
        if self.option('gene2trans').is_set:
            self.logger.debug('{} = {}'.format('gene2trans', self.option('gene2trans').prop['path']))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '16G'

    def end(self):
        super(AnnotFileAgent, self).end()

class AnnotFileTool(Tool):
    def __init__(self, config):
        super(AnnotFileTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.get_isoform2unigene_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/get_isoform2unigene.py')
        self.i2u = os.path.join(self.work_dir, 'i2u')

    def run(self):
        super(AnnotFileTool, self).run()
        if self.option('gtf').is_set and self.option('g2t2p').is_set:
            self.run_get_i2u()
        elif self.option('gene2trans').is_set:
            os.link(self.option('gene2trans').prop['path'], self.i2u)
        self.set_output()
        self.end()

    def run_get_i2u(self):
        cmd = '{} {}'.format(self.python, self.get_isoform2unigene_py)
        cmd += ' -i {}'.format(self.option('gtf').prop['path'])
        cmd += ' -r {}'.format(self.option('g2t2p').prop['path'])
        cmd += ' -o {}'.format(self.i2u)
        cmd_name = 'run_get_i2u'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
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
        source = self.i2u
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.option('i2u', link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test_ref(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_file_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.annot_file',
            'instant': False,
            'options': {
                'gtf': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'g2t2p': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/Annotation_v2/g2t2p',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()