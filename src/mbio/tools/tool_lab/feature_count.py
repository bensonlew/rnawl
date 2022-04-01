# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class FeatureCountAgent(Agent):
    '''
    last_modify: 2020.06.29
    '''

    def __init__(self, parent):
        super(FeatureCountAgent, self).__init__(parent)
        options = [
            {'name': 'sam', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'gtf', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'allowMultiOverlap', 'type': 'bool', 'default': False},
            {'name': 'countMultiMappingReads', 'type': 'bool', 'default': False},
            {'name': 'fraction', 'type': 'bool', 'default': False},
            {'name': 'isPairedEnd', 'type': 'bool', 'default': True},
            {'name': 'count_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(FeatureCountAgent, self).end()


class FeatureCountTool(Tool):
    def __init__(self, config):
        super(FeatureCountTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'Rscript': 'program/R-3.3.1/bin/Rscript'
        }
        self.script = {
            'feature_counts': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/feature_count.r')
        }
        self.file = {
            'count_table': os.path.join(self.output_dir, "count.txt")
        }

    def run(self):
        super(FeatureCountTool, self).run()
        # self.run_sort()
        self.run_feature_count()
        self.set_output()
        self.end()

    def run_feature_count(self):
        cmd = '{} {} {} {} {} {} {} {} {}'.format(self.program['Rscript'], self.script['feature_counts'], self.option('sam').path,
                                      self.option('gtf').path, self.option('isPairedEnd'), self.option('countMultiMappingReads'),
                                                  self.option('allowMultiOverlap'), self.option('fraction'), self.file['count_table'])
        cmd_name = 'run_feature_count'
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
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")


    def set_output(self):
        self.option('count_table').set_path(self.file['count_table'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'feature_count_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.feature_count',
            'instant': False,
            'options': {
                'sam': '/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/HTSeq-0.6.1/scripts/file.sam',
                'gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf',



            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


