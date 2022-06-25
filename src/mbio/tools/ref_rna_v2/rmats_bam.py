# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang,shicaiping,qinjincheng'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import shutil
from mbio.packages.ref_rna_v2.rmats_process_func import process_single_rmats_output_dir
import glob
import unittest

class RmatsBamAgent(Agent):
    '''
    last_modify: 2019.06.13
    '''
    def __init__(self, parent):
        super(RmatsBamAgent, self).__init__(parent)
        options = [
            {'name': 'A_group_bam', 'type': 'string', 'default': None},
            {'name': 'B_group_bam', 'type': 'string', 'default': None},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},
            {'name': 'seq_type', 'type': 'string', 'default': 'paired'},
            {'name': 'lib_type', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'read_length', 'type': 'int', 'default': 150},
            {'name': 'thread', 'type': 'int', 'default': 16},
            {'name': 'tstat', 'type': 'int', 'default': 12},
            {'name': 'cut_off', 'type': 'float', 'default': 0.0001},
        ]
        self.add_option(options)
        self.step.add_steps('rmats_bam')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def step_start(self):
        self.step.rmats_bam.start()
        self.step.update()
    
    def step_end(self):
        self.step.rmats_bam.finish()
        self.step.update()
    
    def set_resource(self):
        self._cpu = self.option('thread')
        self._memory = '80G'

    @toolfuncdeco
    def end(self):
        super(RmatsBamAgent, self).end()

class RmatsBamTool(Tool):
    def __init__(self, config):
        super(RmatsBamTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'rmats': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py')
        }
        self.dir = {
            'output': os.path.join(self.work_dir, 'output_rmats')
        }
        self.file = {
            'ctrl_bam_config': os.path.join(self.dir['output'], 'A_group_bam.txt'),
            'test_bam_config': os.path.join(self.dir['output'], 'B_group_bam.txt')
        }

    @toolfuncdeco
    def run(self):
        super(RmatsBamTool, self).run()
        self.pre_rmats()
        self.run_rmats()
        self.process_output()
        self.set_output()
        self.end()

    @toolfuncdeco
    def pre_rmats(self):
        if os.path.exists(self.dir['output']):
            shutil.rmtree(self.dir['output'])
        os.mkdir(self.dir['output'])
        open(self.file['ctrl_bam_config'], 'w').writelines(self.option('A_group_bam'))
        open(self.file['test_bam_config'], 'w').writelines(self.option('B_group_bam'))

    @toolfuncdeco
    def run_rmats(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats'])
        cmd += ' --gtf {}'.format(self.option('ref_gtf').path)
        cmd += ' --b1 {}'.format(self.file['test_bam_config'])
        cmd += ' --b2 {}'.format(self.file['ctrl_bam_config'])
        cmd += ' --od {}'.format(self.dir['output'])
        cmd += ' -t {}'.format(self.option('seq_type'))
        cmd += ' --libType {}'.format(self.option('lib_type'))
        cmd += ' --readLength {}'.format(self.option('read_length'))
        cmd += ' --nthread {}'.format(self.option('thread'))
        cmd += ' --tstat {}'.format(self.option('tstat'))
        cmd += ' --cstat {}'.format(self.option('cut_off'))
        cmd_name = 'run_rmats'
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
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        elif command.return_code in [-11]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    @toolfuncdeco
    def process_output(self):
        process_single_rmats_output_dir(
            root=self.dir['output'],
            input_gtf=self.option('ref_gtf').path,
            pvalue_fdr='fdr',
            fdr=0.05,
            psi=0
        )

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
    'A5SS.MATS.JC.alter_id.psi_info.txt',
    'A5SS.MATS.JC.alter_id.txt',
    'A5SS.MATS.JCEC.alter_id.psi_info.txt',
    'A5SS.MATS.JCEC.alter_id.txt',
    'A_group_bam.txt',
    'B_group_bam.txt',
    'MXE.MATS.JC.alter_id.psi_info.txt',
    'MXE.MATS.JC.alter_id.txt',
    'MXE.MATS.JCEC.alter_id.psi_info.txt',
    'MXE.MATS.JCEC.alter_id.txt',
    'RI.MATS.JC.alter_id.psi_info.txt',
    'RI.MATS.JC.alter_id.txt',
    'RI.MATS.JCEC.alter_id.psi_info.txt',
    'RI.MATS.JCEC.alter_id.txt',
    'SE.MATS.JC.alter_id.psi_info.txt',
    'SE.MATS.JC.alter_id.txt',
    'SE.MATS.JCEC.alter_id.psi_info.txt',
    'SE.MATS.JCEC.alter_id.txt',
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
    'psi_stats.file.txt'
)


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        A_group_bams = ['/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Con1.bam',
                        '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Con2.bam',
                        '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Con3.bam',
                        '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Con4.bam',
                        '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Con5.bam']
        B_group_bams = ['/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Vit1.bam',
                        '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Vit2.bam',
                        '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Vit3.bam',
                        '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Vit4.bam',
                        '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/bam/Vit5.bam']

        data = {
            'id': 'rmats_bam_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.rmats_bam',
            'instant': False,
            'options': {
                'A_group_bam': ','.join(A_group_bams),
                'B_group_bam': ','.join(B_group_bams),
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/ref_and_new.gtf',
                'seq_type': 'paired',
                'lib_type': 'fr-unstranded',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
