# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang,shicaiping,qinjincheng'

from biocluster.agent import Agent
from mbio.packages.whole_transcriptome.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
from mbio.packages.whole_transcriptome.structure.rmats_process_func import process_single_rmats_output_dir
import glob
import shutil
import unittest

class RmatsAgent(Agent):
    '''
    last_modify: 2019.06.11
    '''
    def __init__(self, parent):
        super(RmatsAgent, self).__init__(parent)
        options = [
            {'name': 'input_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'b1_bam_conf', 'type': 'infile', 'format': 'whole_transcriptome.common'}, # test
            {'name': 'b2_bam_conf', 'type': 'infile', 'format': 'whole_transcriptome.common'}, # ctrl
            {'name': 'read_type', 'type': 'string', 'default': 'paired'},
            {'name': 'lib_type', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'read_length', 'type': 'int', 'default': 150},
            {'name': 'nthread', 'type': 'int', 'default': 10},
            {'name': 'tstat', 'type': 'int', 'default': 10},
            {'name': 'cstat', 'type': 'float', 'default': 0.0001},
        ]
        self.add_option(options)
        self.step.add_steps('rmats')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def step_start(self):
        self.step.rmats.start()
        self.step.update()

    def step_end(self):
        self.step.rmats.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        self._memory = '40G'

    @toolfuncdeco
    def end(self):
        super(RmatsAgent, self).end()

class RmatsTool(Tool):
    def __init__(self, config):
        super(RmatsTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'rmats': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py'),
            'rmats_process': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/structure/rmats_process_script.py')
        }
        self.dir = {
            'output': os.path.join(self.work_dir, 'rmats_output')
        }

    def run(self):
        super(RmatsTool, self).run()
        self.run_rmats()
        self.process_output()
        self.run_rmats_process_script()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_rmats(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats'])
        cmd += ' --gtf {}'.format(self.option('input_gtf').path)
        cmd += ' --b1 {}'.format(self.option('b1_bam_conf').path)
        cmd += ' --b2 {}'.format(self.option('b2_bam_conf').path)
        cmd += ' --od {}'.format(self.dir['output'])
        cmd += ' -t {}'.format(self.option('read_type'))
        cmd += ' --libType {}'.format(self.option('lib_type'))
        cmd += ' --readLength {}'.format(self.option('read_length'))
        cmd += ' --nthread {}'.format(self.option('nthread'))
        cmd += ' --tstat {}'.format(self.option('tstat'))
        cmd += ' --cstat {}'.format(self.option('cstat'))
        cmd_name = 'run_rmats'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def process_output(self):
        process_single_rmats_output_dir(
            root=self.dir['output'],
            input_gtf=self.option('input_gtf').path
        )

    @toolfuncdeco
    def run_rmats_process_script(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats_process'])
        cmd += ' --gtf {}'.format(self.option('input_gtf').path)
        cmd += ' --od {}'.format(self.dir['output'])
        cmd += ' --b1 {}'.format(self.option('b1_bam_conf').path)
        cmd += ' --b2 {}'.format(self.option('b2_bam_conf').path)
        cmd_name = 'run_rmats_process_script'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def set_output(self):
        for p in glob.glob(os.path.join(self.dir['output'], '*')):
            basename = os.path.basename(p)
            if basename in FILES:
                shutil.copy(p, os.path.join(self.output_dir, basename))

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
            'name': 'ref_rna_v3.structure.rmats',
            'instant': False,
            'options': {
                'input_gtf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/ref_and_new.gtf',
                'b1_bam_conf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/b1.bam.conf',
                'b2_bam_conf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/b2.bam.conf',
                'read_type': 'paired',
                'lib_type': 'fr-unstranded'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
