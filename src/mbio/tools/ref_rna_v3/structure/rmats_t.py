# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang,shicaiping,qinjincheng'

import glob
import os
import shutil
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.ref_rna_v3.functions import toolfuncdeco, runcmd
from mbio.packages.ref_rna_v3.structure.rmats_process_func import process_single_rmats_output_dir


class RmatsTAgent(Agent):
    '''
    last_modify: 2019.11.01
    '''

    def __init__(self, parent):
        super(RmatsTAgent, self).__init__(parent)
        options = [
            {'name': 'input_gtf', 'type': 'infile', 'format': 'ref_rna_v3.gtf'},
            {'name': 'b1_bam_conf', 'type': 'infile', 'format': 'ref_rna_v3.common'},  # test
            {'name': 'b2_bam_conf', 'type': 'infile', 'format': 'ref_rna_v3.common'},  # ctrl
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
        self._memory_increase_step = 20

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
        self._cpu = 10
        self._memory = '10G'

    @toolfuncdeco
    def end(self):
        super(RmatsTAgent, self).end()


class RmatsTTool(Tool):
    def __init__(self, config):
        super(RmatsTTool, self).__init__(config)
        self.program = {
                # 'python': 'miniconda2/bin/python',
                'python': 'bioinfo/ref_rna_v3/rmats/miniconda3/bin/python3.7'
        }
        self.script = {
            # 'rmats': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py'),
            # 'rmats': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/medical_transcriptome/rmats-turbo/rmats.py'),
            'rmats': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v3/rmats/miniconda3/rMATS/rmats.py'),
            'rmats_process': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/structure/rmats_process_script.py')
        }
        self.dir = {
            'output': os.path.join(self.work_dir, 'rmats_output'),
            'tmp': os.path.join(self.work_dir, 'tmp_dir')
        }
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.4.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.4.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        super(RmatsTTool, self).run()
        if not self.chk_rmats():
            self.run_rmats()
        self.process_output()
        self.run_rmats_process_script()
        self.set_output()
        self.end()

    def chk_rmats(self):
        if os.path.isdir(self.dir['output']):
            return RESULTS == set(os.listdir(self.dir['output']))
        else:
            return False

    @toolfuncdeco
    def run_rmats(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats'])
        cmd += ' --gtf {}'.format(self.option('input_gtf').path)
        cmd += ' --b1 {}'.format(self.option('b1_bam_conf').path)
        cmd += ' --b2 {}'.format(self.option('b2_bam_conf').path)
        cmd += ' --od {}'.format(self.dir['output'])
        cmd += ' --tmp {}'.format(self.dir['tmp'])
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
            'name': 'ref_rna_v3.structure.rmats_t',
            'instant': False,
            'options': {
                'input_gtf': '/mnt/ilustre/users/isanger/workspace/20201218/Refrna_majorbio_307371/RefrnaAssemble/output/NewTranscripts/ref_and_new.gtf',
                'b1_bam_conf': '/mnt/ilustre/users/isanger/workspace/20201218/Refrna_majorbio_307371/Rmats/Rmats/deltaN.bam.conf',
                'b2_bam_conf': '/mnt/ilustre/users/isanger/workspace/20201218/Refrna_majorbio_307371/Rmats/Rmats/deltaN_IR.bam.conf',
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
