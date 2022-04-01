# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import re
import unittest

class RmatsAgent(Agent):
    '''
    last_modify: 2019.02.22
    '''
    def __init__(self, parent):
        super(RmatsAgent, self).__init__(parent)
        options = [
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            # BAM configuration file of test group
            {'name': 'B1_config', 'type': 'infile', 'format': 'lnc_rna.common'},
            # BAM configuration file of ctrl group
            {'name': 'B2_config', 'type': 'infile', 'format': 'lnc_rna.common'},
            # ['paired', 'single']
            {'name': 'read_type', 'type': 'string', 'default': ''},
            # ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']
            {'name': 'lib_type', 'type': 'string', 'default': ''},
            {'name': 'read_length', 'type': 'int', 'default': 150},
            {'name': 'nthread', 'type': 'int', 'default': 16},
            {'name': 'tstat', 'type': 'int', 'default': 8},
            {'name': 'cstat', 'type': 'float', 'default': 0.0001},
            {'name': 'method', 'type': 'string', 'default': 'fdr'},
            {'name': 'cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'psi', 'type': 'float', 'default': 0.0},
        ]

        self.add_option(options)
        self.step.add_steps('rmats')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 16

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('ref_gtf').is_set:
            self.logger.debug('{} = {}'.format('ref_gtf', self.option('ref_gtf').prop['path']))
            self.infile_size = os.path.getsize(self.option('ref_gtf').prop['path'])
        else:
            raise OptionError('input GTF file must be provided')
        if self.option('B1_config').is_set:
            self.logger.debug('{} = {}'.format('B1_config', self.option('B1_config').prop['path']))
        else:
            raise OptionError('BAM configuration file of test group must be provided')
        if self.option('B2_config').is_set:
            self.logger.debug('{} = {}'.format('B2_config', self.option('B2_config').prop['path']))
        else:
            raise OptionError('BAM configuration file of ctrl group must be provided')
        self.logger.debug('{} = {}'.format('read_type', self.option('read_type')))
        self.logger.debug('{} = {}'.format('lib_type', self.option('lib_type')))
        self.logger.debug('{} = {}'.format('read_length', self.option('read_length')))
        self.logger.debug('{} = {}'.format('nthread', self.option('nthread')))
        self.logger.debug('{} = {}'.format('tstat', self.option('tstat')))
        self.logger.debug('{} = {}'.format('cstat', self.option('cstat')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def step_start(self):
        self.step.rmats.start()
        self.step.update()

    def step_end(self):
        self.step.rmats.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = self.option('nthread')
        self._memory = '40G'

    def end(self):
        super(RmatsAgent, self).end()

class RmatsTool(Tool):
    def __init__(self, config):
        super(RmatsTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.python = 'program/Python/bin/python'
        self.rmats_py = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/rMATS.4.0.2/rMATS-turbo-Linux-UCS2/rmats.py')
        self.process_rmats_output_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/process_rmats_output.py')
        self.rmats_output = os.path.join(self.work_dir, 'rmats_output')
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/rMATS.4.0.1/rMATS-turbo-Linux-UCS2/rmats.py"
        self.dir = {
            'output': os.path.join(self.work_dir, 'rmats_output')
        }

    def chk_rmats(self):
        if os.path.isdir(self.dir['output']):
            return RESULTS == set(os.listdir(self.dir['output']))
        else:
            return False

    def run(self):
        super(RmatsTool, self).run()
        if not self.chk_rmats():
            self.run_rmats()
        self.run_process_rmats_output()
        self.set_output()
        self.end()

    def run_rmats(self):
        with open(self.work_dir + '/rmats.environ', 'w') as f:
            f.write("{}".format(os.environ))

        cmd = '{} {}'.format(self.python, self.rmats_py)
        cmd += ' --gtf {}'.format(self.option('ref_gtf').prop['path'])
        cmd += ' --b1 {}'.format(self.option('B1_config').prop['path'])
        cmd += ' --b2 {}'.format(self.option('B2_config').prop['path'])
        cmd += ' --od {}'.format(self.rmats_output)
        cmd += ' -t {}'.format(self.option('read_type'))
        cmd += ' --libType {}'.format(self.option('lib_type'))
        cmd += ' --readLength {}'.format(self.option('read_length'))
        cmd += ' --nthread {}'.format(self.option('nthread'))
        cmd += ' --tstat {}'.format(self.option('tstat'))
        cmd += ' --cstat {}'.format(self.option('cstat'))
        cmd_name = 'run_rmats'
        self.run_code(cmd_name, cmd)

    def run_process_rmats_output(self):
        cmd = '{} {}'.format(self.python, self.process_rmats_output_py)
        cmd += ' -i {}'.format(self.rmats_output)
        cmd += ' -g {}'.format(self.option('ref_gtf').prop['path'])
        cmd += ' -m {}'.format(self.option('method'))
        cmd += ' -c {}'.format(self.option('cutoff'))
        cmd += ' -p {}'.format(self.option('psi'))
        cmd_name = 'run_process_rmats_output'
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
        result_files = (
            "all_events_detail_big_table.txt",
            "psi_stats.file.txt",
            "event_stats.file.txt",
            "event_type.file.txt"
        )
        for f in os.listdir(self.rmats_output):
            f_path = os.path.join(self.rmats_output, f)
            target = os.path.join(self.output_dir, f)
            if os.path.exists(target):
                os.remove(target)
            if os.path.split(f)[1] in result_files:
                os.link(f_path, target)
            if re.search(r'^fromGTF\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt$', os.path.split(f)[1]):
                os.link(f_path, target)
            if re.search(r'^fromGTF\.novelEvents\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt$', os.path.split(f)[1]):
                os.link(f_path, target)
            if re.search(r'^(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.psi_info\.txt$', os.path.split(f)[1]):
                os.link(f_path, target)
            if re.search(r'^(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.psi_info\.txt$', os.path.split(f)[1]):
                os.link(f_path, target)
            if re.search(r'^(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.txt$', os.path.split(f)[1]):
                os.link(f_path, target)
            if re.search(r'^(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.txt$', os.path.split(f)[1]):
                os.link(f_path, target)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))


RESULTS = {
    'A3SS.MATS.JC.txt',
    'A3SS.MATS.JCEC.txt',
    'A5SS.MATS.JC.txt',
    'A5SS.MATS.JCEC.txt',
    'JC.raw.input.A3SS.txt',
    'JC.raw.input.A5SS.txt',
    'JC.raw.input.MXE.txt',
    'JC.raw.input.RI.txt',
    'JC.raw.input.SE.txt',
    'JCEC.raw.input.A3SS.txt',
    'JCEC.raw.input.A5SS.txt',
    'JCEC.raw.input.MXE.txt',
    'JCEC.raw.input.RI.txt',
    'JCEC.raw.input.SE.txt',
    'MXE.MATS.JC.txt',
    'MXE.MATS.JCEC.txt',
    'RI.MATS.JC.txt',
    'RI.MATS.JCEC.txt',
    'SE.MATS.JC.txt',
    'SE.MATS.JCEC.txt',
    'fromGTF.A3SS.txt',
    'fromGTF.A5SS.txt',
    'fromGTF.MXE.txt',
    'fromGTF.RI.txt',
    'fromGTF.SE.txt',
    'fromGTF.novelEvents.A3SS.txt',
    'fromGTF.novelEvents.A5SS.txt',
    'fromGTF.novelEvents.MXE.txt',
    'fromGTF.novelEvents.RI.txt',
    'fromGTF.novelEvents.SE.txt'
}

FILES = (
    'A3SS.MATS.JC.alter_id.psi_info.txt',
    'A3SS.MATS.JC.alter_id.txt',
    'A3SS.MATS.JCEC.alter_id.psi_info.txt',
    'A3SS.MATS.JCEC.alter_id.txt',
    'A3SS.transcripts.txt',
    'A5SS.MATS.JC.alter_id.psi_info.txt',
    'A5SS.MATS.JC.alter_id.txt',
    'A5SS.MATS.JCEC.alter_id.psi_info.txt',
    'A5SS.MATS.JCEC.alter_id.txt',
    'A5SS.transcripts.txt',
    'MXE.MATS.JC.alter_id.psi_info.txt',
    'MXE.MATS.JC.alter_id.txt',
    'MXE.MATS.JCEC.alter_id.psi_info.txt',
    'MXE.MATS.JCEC.alter_id.txt',
    'MXE.transcripts.txt',
    'RI.MATS.JC.alter_id.psi_info.txt',
    'RI.MATS.JC.alter_id.txt',
    'RI.MATS.JCEC.alter_id.psi_info.txt',
    'RI.MATS.JCEC.alter_id.txt',
    'RI.transcripts.txt',
    'SE.MATS.JC.alter_id.psi_info.txt',
    'SE.MATS.JC.alter_id.txt',
    'SE.MATS.JCEC.alter_id.psi_info.txt',
    'SE.MATS.JCEC.alter_id.txt',
    'SE.transcripts.txt',
    'all_events_detail_big_table.txt',
    'event_stats.file.txt',
    'event_type.file.txt',
    'fromGTF.A3SS.alter_id.txt',
    'fromGTF.A5SS.alter_id.txt',
    'fromGTF.MXE.alter_id.txt',
    'fromGTF.RI.alter_id.txt',
    'fromGTF.SE.alter_id.txt',
    'fromGTF.novelEvents.A3SS.alter_id.txt',
    'fromGTF.novelEvents.A5SS.alter_id.txt',
    'fromGTF.novelEvents.MXE.alter_id.txt',
    'fromGTF.novelEvents.RI.alter_id.txt',
    'fromGTF.novelEvents.SE.alter_id.txt',
    'psi_stats.file.txt',
    'sample.event.coordinates.JCEC.pk',
    'sample.event.coordinates.JC.pk'
)


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.rmats',
            'instant': False,
            'options': {
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'B1_config': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/B1.config',
                'B2_config': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/B2.config',
                'read_type': 'paired',
                'lib_type': 'fr-unstranded',
                'read_length': 150,
                'nthread': 8,
                'tstat': 4,
                'cstat': 0.0001
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
