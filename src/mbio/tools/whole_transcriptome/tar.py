# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
import shutil
import tarfile
import unittest
import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config


class TarAgent(Agent):
    """
    last_modify: 2019.11.01
    """

    def __init__(self, parent):
        super(TarAgent, self).__init__(parent)
        options = [
            {'name': 'indir', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        super(TarAgent, self).end()


class TarTool(Tool):
    def __init__(self, config):
        super(TarTool, self).__init__(config)
        self.pigz = self.config.SOFTWARE_DIR + "/bioinfo/rna/pigz/pigz-2.4/pigz"

    def run(self):
        super(TarTool, self).run()
        self.pack_and_compress()
        self.end()

    def pack_and_compress(self):
        self.logger.warn('start packing folders and compressing them, it will take a long time')
        indir = self.option('indir').prop['path']
        outdir = self.option('indir').prop['path'].replace("output", 'upload')
        self.logger.info("indir: {}".format(indir))
        self.logger.info("outdir: {}".format(outdir))
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        if os.path.exists(indir):
            os.chdir(indir)
        else:
            self.set_error("{}路径不存在".format(indir))
        dirs = os.walk(indir).next()[1]
        for dirname in dirs:
            if os.path.isfile(dirname):
                continue
            tarname = '{}.tar'.format(dirname)
            if dirname in ('file_check', 'rnaseq_mapping', 'gzfastq2fastq'):
                continue
            if not os.path.isfile(os.path.join(indir, tarname)):
                with tarfile.open(tarname, 'w') as tar:
                    tar.add(dirname)
            cmd = self.pigz + ' -p 10 -c ' + tarname + ' > ' + tarname + '.gz'
            self.logger.info(cmd)
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError:
                print("{}压缩出错".format(tarname))
            shutil.rmtree(dirname)
            os.remove(tarname)
        shutil.copytree(indir, outdir)
        os.chdir(self.work_dir)
        self.logger.warn('succeed in packing folders and compressing them')


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'tar_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.tar',
            'instant': False,
            'options': {
                'indir': '/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210323/Longrna_wpm_249763/output/large_gush',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
