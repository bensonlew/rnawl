# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import shutil
import unittest

from biocluster.module import Module


class Macs2Module(Module):
    def __init__(self, work_id):
        super(Macs2Module, self).__init__(work_id)
        options = [
            {'name': 'bam_dir', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'gsize', 'type': 'string', 'default': 'hs'},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} => {}'.format(k, v.value))

    def run(self):
        super(Macs2Module, self).run()
        self.run_macs2()

    def run_macs2(self):
        glob.glob(os.path.join(self.option('bam_dir').path, '*.bam'))
        for bam in glob.glob(os.path.join(self.option('bam_dir').path, '*.bam')):
            macs2 = self.add_tool('whole_transcriptome.macs2')
            treatment = bam
            name = os.path.basename(bam)[:-4]
            opts = {'treatment': bam, 'gsize': self.option('gsize'), 'name': name}
            macs2.set_options(opts)
            self.tools.append(macs2)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for tool in self.tools:
            src = tool.output_dir
            dst = os.path.join(self.output_dir, tool.option('name'))
            shutil.copytree(src, dst)
        self.end()

    def end(self):
        super(Macs2Module, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'macs2_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.macs2',
            'instant': False,
            'options': {
                'bam_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/medip/mapping_data',
                'gsize': 'hs',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
