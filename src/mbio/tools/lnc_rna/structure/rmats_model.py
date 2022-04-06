# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import glob
import unittest

class RmatsModelAgent(Agent):
    '''
    last_modify: 2019.03.18
    '''
    def __init__(self, parent):
        super(RmatsModelAgent, self).__init__(parent)
        options = [
            {'name': 'gene_id', 'type': 'string', 'default': ''},
            {'name': 'event_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'event_type', 'type': 'string', 'default': ''},
            {'name': 'l1', 'type': 'string', 'default': ''}, # test
            {'name': 'l2', 'type': 'string', 'default': ''}, # ctrl
            {'name': 'b1', 'type': 'string', 'default': ''},  # test
            {'name': 'b2', 'type': 'string', 'default': ''},  # ctrl
            {'name': 'exon_s', 'type': 'int', 'default': 1},
            {'name': 'intron_s', 'type': 'int', 'default': 1},
            {'name': 'background', 'type': 'string', 'default': 'white'},
            {'name': 'density', 'type': 'int', 'default': 600},
            {'name': 'quality', 'type': 'int', 'default': 100},
        ]
        self.add_option(options)
        self.step.add_steps('rmats_model')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 16

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def step_start(self):
        self.step.rmats_model.start()
        self.step.update()

    def step_end(self):
        self.step.rmats_model.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        self._memory = '8G'

    def end(self):
        super(RmatsModelAgent, self).end()

class RmatsModelTool(Tool):
    def __init__(self, config):
        super(RmatsModelTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.convert = 'program/ImageMagick/bin/convert'
        self.process_event_file_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/process_event_file.py')
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.rmats2sashimiplot_py = os.path.join(
            self.config.SOFTWARE_DIR, 'bioinfo/rna/rmats2sashimiplot-master/src/rmats2sashimiplot/rmats2sashimiplot.py'
        )
        self.event_file_txt = os.path.join(self.work_dir, 'event_file.txt')
        self.Sashimi_plot = os.path.join(self.work_dir, 'Sashimi_plot')

    def run(self):
        super(RmatsModelTool, self).run()
        self.run_process_event_file()
        self.run_rmats2sashimiplot()
        self.run_convert()
        self.set_output()
        self.end()

    def run_process_event_file(self):
        cmd = '{} {}'.format(self.python, self.process_event_file_py)
        cmd += ' -i {}'.format(self.option('event_file').path)
        cmd += ' -g {}'.format(self.option('gene_id'))
        cmd += ' -o {}'.format(self.event_file_txt)
        cmd_name = 'run_process_event_file'
        self.run_code(cmd_name, cmd)

    def run_rmats2sashimiplot(self):
        cmd = '{} {}'.format(self.python, self.rmats2sashimiplot_py)
        cmd += ' -t {}'.format(self.option('event_type'))
        cmd += ' -e {}'.format(self.event_file_txt)
        cmd += ' --l1 {}'.format(self.option('l1'))
        cmd += ' --l2 {}'.format(self.option('l2'))
        cmd += ' -o {}'.format(self.work_dir)
        cmd += ' --b1 {}'.format(self.option('b1'))
        cmd += ' --b2 {}'.format(self.option('b2'))
        cmd += ' --exon_s {}'.format(self.option('exon_s'))
        cmd += ' --intron_s {}'.format(self.option('intron_s'))
        cmd_name = 'run_rmats2sashimiplot'
        self.run_code(cmd_name, cmd)

    def run_convert(self):
        for n, pdf in enumerate(glob.glob(os.path.join(self.Sashimi_plot, '*.pdf'))):
            cmd = '{}'.format(self.convert)
            cmd += ' -background {}'.format(self.option('background'))
            cmd += ' -density {}'.format(self.option('density'))
            cmd += ' -quality {}'.format(self.option('quality'))
            cmd += ' -flatten {} {}'.format(pdf, '{}.png'.format(pdf[:-4]))
            cmd_name = 'run_convert_{}'.format(n)
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
        for basename in os.listdir(self.Sashimi_plot):
            source = os.path.join(self.Sashimi_plot, basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_list = ['Vit1', 'Vit2', 'Vit3', 'Vit4', 'Vit5']
        ctrl_list = ['Con1', 'Con2', 'Con3', 'Con4', 'Con5']
        dirname = '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/rmats_model'
        l1 = ','.join(test_list)
        l2 = ','.join(ctrl_list)
        b1 = ','.join([os.path.join(dirname, '{}.bam'.format(t)) for t in test_list])
        b2 = ','.join([os.path.join(dirname, '{}.bam'.format(c)) for c in ctrl_list])
        data = {
            'id': 'rmats_model_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.rmats_model',
            'instant': False,
            'options': {
                'gene_id': 'ENSG00000111912',
                'event_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/rmats_model/SE.MATS.JCEC.alter_id.txt',
                'event_type': 'SE',
                'l1': l1,
                'l2': l2,
                'b1': b1,
                'b2': b2,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()