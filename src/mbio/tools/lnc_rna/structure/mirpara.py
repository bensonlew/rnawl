# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import shutil
import unittest

class MirparaAgent(Agent):
    '''
    last_modify: 2019.03.26
    '''
    def __init__(self, parent):
        super(MirparaAgent, self).__init__(parent)
        options = [
            {'name': 'file_in', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'species', 'type': 'string', 'default': ''},
            {'name': 'cutoff', 'type': 'float', 'default': 0.8},
            {'name': 'levels', 'type': 'int', 'default': 20},
            {'name': 'threads', 'type': 'int', 'default': 20},
            {'name': 'file_out', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('mirpara')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.mirpara.start()
        self.step.update()

    def stepfinish(self):
        self.step.mirpara.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('threads')
        self._memory = '8G'

    def end(self):
        super(MirparaAgent, self).end()

class MirparaTool(Tool):
    def __init__(self, config):
        super(MirparaTool, self).__init__(config)
        self.set_environ(PATH=':'.join([
            os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/miRPara/required_packages/unafold-3.8/bin'),
            os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/miRPara/required_packages/ct2out'),
            os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/miRPara/required_packages/libsvm-3.14')
        ]))
        self.python = 'miniconda2/bin/python'
        self.perl = 'bioinfo/lnc_rna/miniconda2/bin/perl'
        self.mirpara_pretreat_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/mirpara_pretreat.py')
        self.mirpara_pl = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/miRPara/miRPara.pl')
        self.mirpara_posttest_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/mirpara_posttest.py')
        self.map_txt = os.path.join(self.work_dir, 'map.txt')
        self.file_in = os.path.join(self.work_dir, os.path.basename(self.option('file_in').path))
        self.mirpara_out = '{}_level_{}.out'.format(self.file_in[:-3], self.option('levels'))
        self.file_out = os.path.join(self.work_dir, 'mirpara.tblout')

    def run(self):
        super(MirparaTool, self).run()
        self.run_mirpara_pretreat()
        self.run_mirpara()
        self.run_mirpara_posttest()
        self.set_output()
        self.end()

    def run_mirpara_pretreat(self):
        cmd = '{} {}'.format(self.python, self.mirpara_pretreat_py)
        cmd += ' -i {}'.format(self.option('file_in').path)
        cmd += ' -m {}'.format(self.map_txt)
        cmd += ' -o {}'.format(self.file_in)
        cmd_name = 'run_mirpara_pretreat'
        self.run_code(cmd_name, cmd)

    def run_mirpara(self):
        cmd = '{} {}'.format(self.perl, self.mirpara_pl)
        cmd += ' -s {}'.format(self.option('species'))
        cmd += ' -c {}'.format(self.option('cutoff'))
        cmd += ' -l {}'.format(self.option('levels'))
        cmd += ' -t {}'.format(self.option('threads'))
        cmd += ' {}'.format(self.file_in)
        cmd_name = 'run_mirpara'
        self.run_code(cmd_name, cmd)

    def run_mirpara_posttest(self):
        cmd = '{} {}'.format(self.python, self.mirpara_posttest_py)
        cmd += ' -i {}'.format(self.mirpara_out)
        cmd += ' -m {}'.format(self.map_txt)
        cmd += ' -o {}'.format(self.file_out)
        cmd_name = 'run_mirpara_posttest'
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
        source = self.file_out
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('file_out').set_path(link_name)
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
            'id': 'mirpara_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.mirpara',
            'instant': False,
            'options': {
                'file_in': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/precursor_candidate/output/candidate.fa',
                'species': 'hsa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()