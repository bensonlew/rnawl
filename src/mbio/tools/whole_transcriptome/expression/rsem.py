# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class RsemAgent(Agent):
    '''
    last_modify: 2019.11.18
    '''

    def __init__(self, parent):
        super(RsemAgent, self).__init__(parent)
        STRANDEDNESS = ('none', 'forward', 'reverse')
        options = [
            {'name': 'paired_end', 'type': 'bool', 'default': True},
            {'name': 'upstream_read_file', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'downstream_read_file', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'reference_name', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'strandedness', 'type': 'string', 'default': STRANDEDNESS[2]},
            {'name': 'threads', 'type': 'int', 'default': 8},
            {'name': 'estimate_rspd', 'type': 'bool', 'default': True},
            {'name': 'isoforms_results', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'genes_results', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '50G'

    def end(self):
        super(RsemAgent, self).end()


class RsemTool(Tool):
    def __init__(self, config):
        super(RsemTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.program = {
            'rsem_calculate_expression': 'bioinfo/rna/RSEM-1.3.1/rsem-calculate-expression'
        }
        self.file = {
            'isoforms_results': os.path.join(self.work_dir, '{}.isoforms.results'.format(self.option('sample'))),
            'genes_results': os.path.join(self.work_dir, '{}.genes.results'.format(self.option('sample')))
        }
        self.dir = {
            'bowtie2_path': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/miniconda3/bin'),
            'glibc_2.14': os.path.join(self.config.SOFTWARE_DIR, 'library/glibc-2.14/lib')
        }
        self.set_environ(PATH=self.dir['bowtie2_path'])
        self.set_environ(LD_LIBRARY_PATH=self.dir['glibc_2.14'])

    def run(self):
        super(RsemTool, self).run()
        self.pre_rsem()
        self.run_rsem_calculate_expression()
        self.set_output()
        self.end()

    def pre_rsem(self):
        import subprocess
        proc = subprocess.Popen('which samtools', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc.wait()
        self.logger.info(proc.communicate())
        proc = subprocess.Popen('echo ${LD_LIBRARY_PATH}', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc.wait()
        self.logger.info(proc.communicate())
        proc = subprocess.Popen('echo ${PATH}', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc.wait()
        self.logger.info(proc.communicate())

    def run_rsem_calculate_expression(self):
        if self.option("downstream_read_file").is_set:
            self.option("paired_end", True)
        else:
            self.option("paired_end", False)
        self.logger.info("paired_end: {}".format(self.option("paired_end")))
        cmd = '{}'.format(self.program['rsem_calculate_expression'])
        cmd += ' --bowtie2 --bowtie2-path {}'.format(self.dir['bowtie2_path'])
        if self.option('estimate_rspd'):
            cmd += ' --estimate-rspd'
        cmd += ' --num-threads {}'.format(self.option('threads'))
        cmd += ' --strandedness {}'.format(self.option('strandedness'))
        if self.option('paired_end'):
            cmd += ' --paired-end {} {} {} {}'.format(
                self.option('upstream_read_file').path,
                self.option('downstream_read_file').path,
                self.option('reference_name').path,
                self.option('sample'))
        else:
            cmd += ' {} {} {}'.format(
                self.option('upstream_read_file').path,
                self.option('reference_name').path,
                self.option('sample'))
        runcmd(self, 'run_rsem_calculate_expression', cmd)

    def set_output(self):
        for name in ('isoforms_results', 'genes_results'):
            source = self.file[name]
            link_name = os.path.join(self.output_dir, os.path.basename(source))
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.option(name).set_path(link_name)


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rsem_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.rsem',
            'instant': False,
            'options': {
                'upstream_read_file': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/FastpRna/output/fastq/C_6w.clean.1.fastq',
                'downstream_read_file': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/FastpRna/output/fastq/C_6w.clean.2.fastq',
                'reference_name': '/mnt/ilustre/users/sanger-dev/workspace/20191118/Single_index_6600_3886/Index/output/rsem.index',
                'sample': 'C_6w',
                'strandedness': 'reverse'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
