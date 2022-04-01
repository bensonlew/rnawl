# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang,shicaiping,qinjincheng'

from biocluster.module import Module
from mbio.packages.whole_transcriptome.functions import modlfuncdeco
import os
import glob
import shutil
import unittest

class RmatsModule(Module):
    '''
    last_modify: 2019.06.12
    '''
    def __init__(self, work_id):
        super(RmatsModule, self).__init__(work_id)
        options = [
            {'name': 'control_table', 'type': 'infile', 'format': 'sample.control_table'},
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'bam_input', 'type': 'infile', 'format': 'align.bwa.bam_dir'},
            {'name': 'input_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            # ['paired', 'single']
            {'name': 'seq_type', 'type': 'string', 'default': 'paired'},
            # ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']
            {'name': 'lib_type', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'read_length', 'type': 'int', 'default': 150},
            {'name': 'cstat', 'type': 'float', 'default': 0.0001},
        ]
        self.add_option(options)

    @modlfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @modlfuncdeco
    def run(self):
        super(RmatsModule, self).run()
        self.run_rmats()

    @modlfuncdeco
    def run_rmats(self):
        num, vs_list = self.option('control_table').get_control_info()
        self.group2samples = self.option('group_table').get_group_spname()
        self.sample2path = {os.path.basename(p).rstrip('.bam'): p
                            for p in glob.glob(os.path.join(self.option('bam_input').path, '*.bam'))}
        self.tools = [self.set_rmats(n, ctrl, test) for n, (ctrl, test) in enumerate(vs_list)]
        if len(self.tools) == 1:
            self.tools[0].on('end', self.run_rmats_count)
        else:
            self.on_rely(self.tools, self.run_rmats_count)
        for tool in self.tools:
            tool.run()

    @modlfuncdeco
    def set_rmats(self, n, ctrl, test):
        self.step.add_steps('rmats_{}'.format(n))
        rmats = self.add_tool('whole_transcriptome.structure.rmats')
        b1_bam_conf = os.path.join(rmats.work_dir, '{}.bam.conf'.format(test))
        b2_bam_conf = os.path.join(rmats.work_dir, '{}.bam.conf'.format(ctrl))
        open(b1_bam_conf, 'w').write(','.join(self.sample2path[s] for s in self.group2samples[test]))
        open(b2_bam_conf, 'w').write(','.join(self.sample2path[s] for s in self.group2samples[ctrl]))
        rmats.set_options({
            'input_gtf': self.option('input_gtf').path,
            'b1_bam_conf': b1_bam_conf,
            'b2_bam_conf': b2_bam_conf,
            'read_type': self.option('seq_type'),
            'lib_type': self.option('lib_type'),
            'read_length': self.option('read_length'),
            'cstat': self.option('cstat')
        })
        rmats.vs_name = '{}_vs_{}'.format(ctrl, test)
        rmats.on('start', self.set_step, {'start': getattr(self.step, 'rmats_{}'.format(n))})
        rmats.on('end', self.set_step, {'end': getattr(self.step, 'rmats_{}'.format(n))})
        return rmats

    def run_rmats_count(self):
        self.step.add_steps('rmats_count')
        self.rmats_count = self.add_tool('whole_transcriptome.structure.rmats_count')
        jcpklist = os.path.join(self.rmats_count.work_dir, 'jc.pk.list')
        jcecpklist = os.path.join(self.rmats_count.work_dir, 'jcec.pk.list')
        open(jcpklist, 'w').writelines(
            '{}\n'.format(os.path.join(tool.output_dir, 'sample.event.coordinates.JC.pk')) for tool in self.tools
        )
        open(jcecpklist, 'w').writelines(
            '{}\n'.format(os.path.join(tool.output_dir, 'sample.event.coordinates.JCEC.pk')) for tool in self.tools
        )
        self.rmats_count.set_options({
            'jcpklist': jcpklist,
            'jcecpklist': jcecpklist
        })
        self.rmats_count.on('start', self.set_step, {'start': getattr(self.step, 'rmats_count')})
        self.rmats_count.on('end', self.set_step, {'end': getattr(self.step, 'rmats_count')})
        self.rmats_count.on('end', self.set_output)
        self.rmats_count.run()

    def set_output(self):
        shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir)
        for tool in self.tools:
            shutil.copytree(tool.output_dir, os.path.join(self.output_dir, tool.vs_name))
        else:
            shutil.copy(self.rmats_count.option('ojc').path,
                        os.path.join(self.output_dir, os.path.basename(self.rmats_count.option('ojc').path)))
            shutil.copy(self.rmats_count.option('ojcec').path,
                        os.path.join(self.output_dir, os.path.basename(self.rmats_count.option('ojcec').path)))
            self.end()

    @modlfuncdeco
    def end(self):
        super(RmatsModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run script to do test.
    '''

    def test(self):
        import random
        from biocluster.wsheet import Sheet
        from mbio.workflows.single import SingleWorkflow
        data = {
            'id': 'rmats_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.rmats',
            'instant': False,
            'options': {
                'control_table': '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/whole_transcriptome/rmats/test_data/test_control.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/whole_transcriptome/rmats/test_data/test_group.txt',
                'bam_input': '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/whole_transcriptome/rmats/test_data/test',
                'input_gtf': '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/whole_transcriptome/rmats/test_data/ref_and_new.gtf',
                'seq_type': 'paired',
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
