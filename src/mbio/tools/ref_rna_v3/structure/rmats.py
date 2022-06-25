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


class RmatsAgent(Agent):
    '''
    last_modify: 2019.11.01
    '''

    def __init__(self, parent):
        super(RmatsAgent, self).__init__(parent)
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
            {'name': 'version', 'type': 'string', 'default': 'rMATS.4.0.2'},
        ]
        self.add_option(options)
        self.step.add_steps('rmats')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 40

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
        self._cpu = self.option("nthread")
        self._memory = '10G'

    @toolfuncdeco
    def end(self):
        super(RmatsAgent, self).end()


class RmatsTool(Tool):
    def __init__(self, config):
        super(RmatsTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        if self.option("version") == "rMATS.4.0.2":
            self.program = {
                'python': 'miniconda2/bin/python',
            }
            self.script = {
                'rmats': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py'),
                'rmats_process': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/structure/rmats_process_script.py')
            }
            self.dir = {
                'output': os.path.join(self.work_dir, 'rmats_output')
            }
        elif self.option("version") == "rMATS.4.1.1":
            self.program = {
                'python3': 'bioinfo/ref_rna_v3/rmats_4.1.1/miniconda3/bin/python',
                'python': 'miniconda2/bin/python',
            }
            self.script = {
                'rmats': os.path.join(self.config.SOFTWARE_DIR,
                                      'bioinfo/ref_rna_v3/rmats_4.1.1/miniconda3/bin/rmats.py'),
                'rmats_process': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/structure/rmats_process_script.py')
            }
            self.dir = {
                'output': os.path.join(self.work_dir, 'rmats_output'),
                'tmp_output': os.path.join(self.work_dir, 'rmats_tmp_output'),
            }

    def run(self):
        super(RmatsTool, self).run()
        if not self.chk_rmats():
            self.run_rmats()
        self.process_output()
        self.run_rmats_process_script()
        self.set_output()
        self.end()

    def chk_rmats(self):
        if os.path.isdir(self.dir['output']):
            if self.option("version") == 'rMATS.4.1.1':
                if RESULTS_new.issubset(set(os.listdir(self.dir['output']))):
                    return True
            else:
                return RESULTS == set(os.listdir(self.dir['output']))
        else:
            return False

    @toolfuncdeco
    def run_rmats(self):
        if os.path.exists(self.dir['output']):
            shutil.rmtree(self.dir['output'])
        if self.option("version") == "rMATS.4.1.1":
            if os.path.exists(self.dir['tmp_output']):
                shutil.rmtree(self.dir['tmp_output'])
            cmd = '{} {}'.format(self.program['python3'], self.script['rmats'])
        else:
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
        if self.option("version") == "rMATS.4.1.1":
            cmd += ' --variable-read-length'
            #cmd += ' --allow-clipping'
            cmd += ' --tmp {}'.format(self.dir['tmp_output'])
            cmd += ' --novelSS'
        cmd_name = 'run_rmats'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("rmats运行完成")
        elif command.return_code in [2,255,-11]:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("rmats运行出错, Tool退出")

    @toolfuncdeco
    def process_output(self):
        if self.option("version") == 'rMATS.4.1.1':
            events = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']
            for event in events:
                novelJunction = os.path.join(self.dir['output'], 'fromGTF.novelJunction.{}.txt'.format(event))
                novelSpliceSite = os.path.join(self.dir['output'], 'fromGTF.novelJunction.{}.txt'.format(event))
                if not (os.path.exists(novelJunction) and os.path.exists(novelSpliceSite)):
                    self.set_error("rmats运行结果文件不完整")
                with open(novelJunction, "r") as f1, open(novelSpliceSite, "r") as f2, open(os.path.join(self.dir['output'], 'fromGTF.novelEvents.{}.txt'.format(event)), "w") as w1:
                    for line in f1:
                        w1.write(line)
                    head = f2.readline()
                    for line in f2:
                        w1.write(line)
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
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        for run_times in range(1,51):
            self.logger.info('运行rmats软件第{}次，进行可变剪接分析'.format(run_times))
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("rmats运行完成")
                break
            elif command.return_code in [2,255]:
                self.logger.info("return code: {}".format(command.return_code))
                self.add_state('memory_limit', 'memory is low!')
            else:
                self.logger.info("第{}次尝试，返回值为{}".format(run_times, command.return_code))
                if run_times == 50:
                    self.set_error("rmats多次运行都没有完成，退出!", code = "32006109")
                else:
                    self.logger.info("rmats运行出错,进行下一次尝试")


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

RESULTS_new = {
    'fromGTF.novelSpliceSite.RI.txt',
    'fromGTF.novelSpliceSite.A3SS.txt',
    'JCEC.raw.input.A3SS.txt',
    'JCEC.raw.input.MXE.txt',
    'A3SS.MATS.JCEC.txt',
    'JCEC.raw.input.A5SS.txt',
    'JC.raw.input.RI.txt',
    'tmp',
    'SE.MATS.JCEC.txt',
    'MXE.MATS.JCEC.txt',
    'JCEC.raw.input.SE.txt',
    'fromGTF.novelSpliceSite.SE.txt',
    'fromGTF.novelJunction.A5SS.txt',
    'SE.MATS.JC.txt',
    'RI.MATS.JCEC.txt',
    'RI.MATS.JC.txt',
    'A5SS.MATS.JCEC.txt',
    'summary.txt',
    'fromGTF.novelJunction.A3SS.txt',
    'A3SS.MATS.JC.txt',
    'fromGTF.novelSpliceSite.MXE.txt',
    'fromGTF.novelJunction.MXE.txt',
    'fromGTF.A3SS.txt',
    'JC.raw.input.A5SS.txt',
    'JC.raw.input.A3SS.txt',
    'JC.raw.input.SE.txt',
    'fromGTF.MXE.txt',
    'JC.raw.input.MXE.txt',
    'fromGTF.RI.txt',
    'fromGTF.novelJunction.SE.txt',
    'fromGTF.A5SS.txt',
    'fromGTF.SE.txt',
    'fromGTF.novelSpliceSite.A5SS.txt',
    'JCEC.raw.input.RI.txt',
    'fromGTF.novelJunction.RI.txt',
    'MXE.MATS.JC.txt',
    'A5SS.MATS.JC.txt'
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
