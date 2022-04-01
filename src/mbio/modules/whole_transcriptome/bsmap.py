# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import shutil
import gevent.subprocess as subprocess
import unittest

from biocluster.module import Module
from mbio.packages.whole_transcriptome.utils import read_fastq_dir


class BsmapModule(Module):
    def __init__(self, work_id):
        super(BsmapModule, self).__init__(work_id)
        options = [
            {'name': 'reference', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'combine', 'type': 'bool', 'default': False},
        ]
        self.add_option(options)
        self.tools = list()
        self.file = {'ratio': os.path.join(self.output_dir, 'meth.txt')}
        self.script = {'methratio': os.path.join(self.get_workflow().config.SOFTWARE_DIR,
                                                 'bioinfo/wgbs/miniconda3/bin/methratio.py')}
        self.program = {'python': os.path.join(self.get_workflow().config.SOFTWARE_DIR, 'program/Python/bin/python')}

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} => {}'.format(k, v.value))

    def run(self):
        super(BsmapModule, self).run()
        self.run_bsmap()

    def run_bsmap(self):
        is_se, fastqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        for sample, fastqs in fastqs_dict.items():
            bsmap = self.add_tool('whole_transcriptome.bsmap')
            opts = {'query_a': fastqs[0], 'query_b': fastqs[1], 'reference': self.option('reference').path,
                    'sample': sample}
            bsmap.set_options(opts)
            self.tools.append(bsmap)
        else:
            if self.option('combine'):
                self.on_rely(self.tools, self.run_methratio)
            else:
                self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def run_methratio(self):
        cmd = '{} {}'.format(self.program['python'], self.script['methratio'])
        cmd += ' -o {}'.format(self.file['ratio'])
        cmd += ' -d {}'.format(self.option('reference').path)
        for tool in self.tools:
            cmd += ' {}'.format(os.path.join(tool.output_dir, '{}.bsp'.format(tool.option('sample'))))
        self.logger.debug(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             env=os.environ, universal_newlines=True, bufsize=0)
        stdoutdata, stderrdata = p.communicate()
        self.logger.debug(stdoutdata)
        self.logger.debug(stderrdata)
        p.stdout.close()
        returncode = p.poll()
        if returncode:
            self.set_error(stderrdata)
        else:
            self.set_output()

    def set_output(self):
        for path in glob.glob(os.path.join(self.output_dir, '*')):
            if os.path.isdir(path):
                shutil.rmtree(path)
        for tool in self.tools:
            src = tool.output_dir
            dst = os.path.join(self.output_dir, tool.option('sample'))
            shutil.copytree(src, dst)
        self.end()

    def end(self):
        super(BsmapModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the module. Just run this script to do test.
    """

    def test_rcrs(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bsmap_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.bsmap',
            'instant': False,
            'options': {
                'reference': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens'
                             '/Mitochondrion/rCRS/rCRS.fasta',
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
            'id': 'bsmap_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.bsmap',
            'instant': False,
            'options': {
                'reference': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens'
                             '/Mitochondrion/RSRS/RSRS.fasta',
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
