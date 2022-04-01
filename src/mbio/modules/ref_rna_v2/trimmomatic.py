# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
import os
import glob
import unittest

class TrimmomaticModule(Module):
    '''
    last_modify: 2019.04.29
    '''
    def __init__(self, work_id):
        super(TrimmomaticModule, self).__init__(work_id)
        options = [
            {'name': 'fq_dir', 'type': 'infile', 'format': 'ref_rna_v2.fastq_dir'},
            {'name': 'fq_type', 'type': 'string', 'default': None}, # PE SE
            {'name': 'phred', 'type': 'int', 'default': None}, # 33 64
            {'name': 'fq_list', 'type': 'outfile', 'format': 'ref_rna_v2.fastq_list'},
        ]
        self.add_option(options)
        self.tools = list()
    
    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.option('fq_dir').prepare()
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(TrimmomaticModule, self).run()
        self.run_trimmomatic()

    def run_trimmomatic(self):
        for n, sample in enumerate(self.option('fq_dir').samples):
            self.step.add_steps('trimmomatic_{}'.format(n))
            trimmomatic = self.add_tool('ref_rna_v2.trimmomatic')
            options = {
                'fq_type': self.option('fq_type'),
                'sample_name': sample,
            }
            if self.option('fq_type') == 'PE':
                options.update({
                    'input_file1': self.option('fq_dir').pe1_fastq[sample],
                    'input_file2': self.option('fq_dir').pe2_fastq[sample],
                })
            elif self.option('fq_type') == 'SE':
                options.update({'input_file': self.option('fq_dir').se_fastq[sample]})
            trimmomatic.set_options(options)
            trimmomatic.on('start', self.set_step, {'start': getattr(self.step, 'trimmomatic_{}'.format(n))})
            trimmomatic.on('end', self.set_step, {'end': getattr(self.step, 'trimmomatic_{}'.format(n))})
            self.tools.append(trimmomatic)
        else:
            if len(self.tools) == 1:
                self.tools[0].on('end', self.set_output)
            else:
                self.on_rely(self.tools, self.set_output)
            for tool in self.tools:
                tool.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        f2d_lines = list()
        n2l_lines = list()
        for tool in self.tools:
            if self.option('fq_type') == 'PE':
                fq1_path = tool.option('output_file_1p').path
                fq2_path = tool.option('output_file_2p').path
                f2d_lines.append('{}\t{}\t{}\n'.format(os.path.basename(fq1_path), tool.option('sample_name'), 'l'))
                f2d_lines.append('{}\t{}\t{}\n'.format(os.path.basename(fq2_path), tool.option('sample_name'), 'r'))
                n2l_lines.append('{}\t{}\t{}\n'.format(tool.option('sample_name'), fq1_path, fq2_path))
            elif self.option('fq_type') == 'SE':
                fq_path = tool.option('output_file').path
                f2d_lines.append('{}\t{}\n'.format(os.path.basename(fq_path), tool.option('sample_name')))
                n2l_lines.append('{}\t{}\n'.format(tool.option('sample_name'), fq_path))
            for source in glob.glob(os.path.join(tool.output_dir, '*')):
                link_name = os.path.join(self.output_dir, os.path.basename(source))
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
                self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        else:
            f2d_file = os.path.join(self.output_dir, 'list.txt')
            n2l_file = os.path.join(self.output_dir, 'fq.list')
            open(f2d_file, 'w').writelines(f2d_lines)
            open(n2l_file, 'w').writelines(n2l_lines)
            self.option('fq_list').set_path(n2l_file)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(TrimmomaticModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_pe(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'trimmomatic_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'ref_rna_v2.trimmomatic',
            'instant': False,
            'options': {
                'fq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_demo/remote_input/fastq_dir',
                'fq_type': 'PE',
                'phred': 33,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_se(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'trimmomatic_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'ref_rna_v2.trimmomatic',
            'instant': False,
            'options': {
                'fq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/chip/raw_data',
                'fq_type': 'SE',
                'phred': 33,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_pe')])
    unittest.TextTestRunner(verbosity=2).run(suite)
