# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class HtseqCountAgent(Agent):
    '''
    last_modify: 2020.06.29
    '''

    def __init__(self, parent):
        super(HtseqCountAgent, self).__init__(parent)
        options = [
            {'name': 'sam', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'bamorsam', 'type': 'string', 'default': 'sam'},
            {'name': 'order', 'type': 'string', 'default': 'name'},
            {'name': 'gtf', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'method', 'type': 'string', 'default':'union'},
            {'name': 'nonunique', 'type': 'string', 'default': 'none'},
            {'name': 'gort', 'type': 'string', 'default': 'gene_id'},
            {'name': 'sample', 'type': 'string'},
            {'name': 'count_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(HtseqCountAgent, self).end()


class HtseqCountTool(Tool):
    def __init__(self, config):
        super(HtseqCountTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'python': 'miniconda2/bin/python',
            # 'python': 'miniconda2/bin/python',
            # 'htseq_count': 'bioinfo/rna/HTSeq-0.6.1/scripts/htseq-count',
            # 'htseq_count': 'miniconda2/bin/htseq-count',
            # 'python': 'bioinfo/ref_rna_v3/HTSeq/miniconda2/bin/python',
            'htseq_count': 'bioinfo/ref_rna_v3/HTSeq/miniconda2/bin/htseq-count',
            'samtools': '/bioinfo/align/samtools-1.3.1/samtools'
        }
        # self.script = {
        #     'table_kit': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/table_kit.py')
        # }
        self.file = {
            'sort_': os.path.join(self.work_dir, 'sort.sam'),
            'count_table': os.path.join(self.output_dir, "{}_count.txt".format(self.option('sample')))
        }

    def run(self):
        super(HtseqCountTool, self).run()
        # self.run_sort()
        self.run_htseq_count()
        self.set_output()
        self.end()

    def run_sort(self):
        cmd = '{} sort -n {} > {}'.format(self.program['samtools'], self.option('sam').path, self.file['sort_'])
        cmd_name = 'run_samtools'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_htseq_count(self):
        cmd = '{} -f {} -r {} '.format(self.program['htseq_count'], self.option('bamorsam'), self.option('order'))
        cmd += '-a 10 -t exon -i {} -m {} --stranded=no --nonunique {} '.format(self.option('gort'), self.option('method'),
                                                                                self.option('nonunique'))
        # cmd += '-a 10 -t exon -i {} -m {} --stranded=no '.format(self.option('gort'), self.option('method'))
        cmd += '{} {} > {}'.format(self.option('sam').path, self.option('gtf').path, self.file['count_table'])
        cmd_name = 'run_htseq_count'
        runcmd(self, cmd_name, cmd, shell=True)

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
            'id': 'htseq_count_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.htseq_count',
            'instant': False,
            'options': {
                'sam': '/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/HTSeq-0.6.1/scripts/file.sam',
                'gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf',
                'method': 'union',
                'sample': 'zjx',
                'nonunique': 'all'


            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


