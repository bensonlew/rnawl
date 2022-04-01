# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.module import Module
from mbio.packages.whole_transcriptome.utils import read_fastq_dir


class BismarkModule(Module):
    def __init__(self, work_id):
        super(BismarkModule, self).__init__(work_id)
        options = [
            {'name': 'genome_folder', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} => {}'.format(k, v.value))

    def run(self):
        super(BismarkModule, self).run()
        self.run_bismark()

    def run_bismark(self):
        is_se, fastqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        for sample, fastqs in fastqs_dict.items():
            bismark = self.add_tool('whole_transcriptome.bismark')
            opts = {'sample': sample, 'genome_folder': self.option('genome_folder').path,
                    'mates1': fastqs[0], 'mates2': fastqs[1]}
            bismark.set_options(opts)
            self.tools.append(bismark)
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
            dst = os.path.join(self.output_dir, tool.option('sample'))
            shutil.copytree(src, dst)
        self.end()

    def end(self):
        super(BismarkModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the module. Just run this script to do test.
    """

    def test_rcrs(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bismark_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.bismark',
            'instant': False,
            'options': {
                'genome_folder': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates'
                                 '/Homo_sapiens/Mitochondrion/rCRS',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/wgbs/clean_data',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_rsrs(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bismark_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.bismark',
            'instant': False,
            'options': {
                'genome_folder': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates'
                                 '/Homo_sapiens/Mitochondrion/rCRS',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/wgbs/clean_data',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_rcrs'), TestFunction('test_rsrs')])
    unittest.TextTestRunner(verbosity=2).run(suite)
