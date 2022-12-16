# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import tarfile
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool


class TransferAgent(Agent):
    """
    last_modify: 2019.11.01
    """

    def __init__(self, parent):
        super(TransferAgent, self).__init__(parent)
        options = [
            {'name': 'indir', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'outdir', 'type': 'outfile', 'format': 'whole_transcriptome.common_dir'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(TransferAgent, self).end()


class TransferTool(Tool):
    def __init__(self, config):
        super(TransferTool, self).__init__(config)

    def run(self):
        super(TransferTool, self).run()
        self.set_output()
        self.decompress()
        self.end()

    def set_output(self):
        self.logger.info(self.option('indir').path)
        for name in os.listdir(self.option('indir').path):
            src = os.path.join(self.option('indir').path, name)
            dst = os.path.join(self.output_dir, name)
            self.logger.info(name)
            self.logger.info(src)
            self.logger.info(dst)
            if os.path.isfile(src):
                self.link2outputdir(src, dst)
            elif os.path.isdir(src):
                self.move2outputdir(src, name)
        self.option('outdir').set_path(self.output_dir)

    def decompress(self):
        self.logger.warn('start decompressing compressed file, it will take a long time')
        os.chdir(self.output_dir)
        nondirs = os.walk(self.output_dir).next()[2]
        for fname in nondirs:
            if fname.endswith('.tar.gz'):
                self.logger.debug('start decompressing {}'.format(fname))
                with tarfile.open(fname, 'r:gz') as tar:
                    def is_within_directory(directory, target):
                        
                        abs_directory = os.path.abspath(directory)
                        abs_target = os.path.abspath(target)
                    
                        prefix = os.path.commonprefix([abs_directory, abs_target])
                        
                        return prefix == abs_directory
                    
                    def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                    
                        for member in tar.getmembers():
                            member_path = os.path.join(path, member.name)
                            if not is_within_directory(path, member_path):
                                raise Exception("Attempted Path Traversal in Tar File")
                    
                        tar.extractall(path, members, numeric_owner=numeric_owner) 
                        
                    
                    safe_extract(tar, self.output_dir)
                self.logger.debug('succeed in decompressing {}'.format(fname))
        os.chdir(self.work_dir)
        self.logger.warn('succeed in decompressing compressed file')

    def link2outputdir(self, src, dst):
        if os.path.isfile(dst):
            os.remove(dst)
        os.link(src, dst)
        self.logger.debug('succeed in linking {} to {}'.format(src, dst))

    def move2outputdir(self, olddir, newname):
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        else:
            self.logger.debug('succeed in linking {} to {}'.format(olddir, newdir))

    def move_file(self, src, dst):
        if os.path.isfile(src):
            os.link(src, dst)
        else:
            os.mkdir(dst)
            for _file in os.listdir(src):
                old_path = os.path.join(src, _file)
                new_path = os.path.join(dst, _file)
                self.move_file(old_path, new_path)


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'transfer_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.transfer',
            'instant': False,
            'options': {
                'indir': 's3nb://refrnav2/files/m_34394/34394_5d80a9c6283a7/i-sanger_204579/intermediate_results'
                         '/Align/AlignBam/',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
