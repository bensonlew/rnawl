# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd


class BismarkAgent(Agent):
    """
    last_modify: 2020.02.18
    """

    def __init__(self, parent):
        super(BismarkAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'nucleotide_coverage', 'type': 'bool', 'default': True},
            {'name': 'genome_folder', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'mates1', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'mates2', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'non_directional', 'type': 'bool', 'default': True},
            {'name': 'comprehensive', 'type': 'bool', 'default': True},
            {'name': 'parallel', 'type': 'int', 'default': 4},
            {'name': 'bed_graph', 'type': 'bool', 'default': True},
            {'name': 'cytosine_report', 'type': 'bool', 'default': True},
            {'name': 'cx_context', 'type': 'bool', 'default': True},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} => {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(BismarkAgent, self).end()


class BismarkTool(Tool):
    def __init__(self, config):
        super(BismarkTool, self).__init__(config)
        self.program = {
            'bismark': 'bioinfo/wgbs/miniconda3/bin/bismark',
            'deduplicate_bismark': 'bioinfo/wgbs/miniconda3/bin/deduplicate_bismark',
            'bismark_methylation_extractor': 'bioinfo/wgbs/miniconda3/bin/bismark_methylation_extractor',
            'bismark2report': 'bioinfo/wgbs/miniconda3/bin/bismark2report',
        }
        self.dir = {
            'result': os.path.join(self.output_dir, 'result'),
            'report': os.path.join(self.output_dir, 'report'),
        }
        self.file = {
            'bam': os.path.join(self.dir['result'], '{}_pe.bam'.format(self.option('sample'))),
            'dedup_bam': os.path.join(self.dir['result'], '{}_pe.deduplicated.bam'.format(self.option('sample')))
        }
        for path in self.dir.values():
            if os.path.isdir(path):
                shutil.rmtree(path)
            os.mkdir(path)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/wgbs/miniconda3/bin'))

    def run(self):
        super(BismarkTool, self).run()
        self.run_bismark()
        self.run_deduplicate_bismark()
        self.run_bismark_methylation_extractor()
        self.run_bismark2report()
        self.set_output()
        self.end()

    def run_bismark(self):
        cmd = '{} -un --basename {}'.format(self.program['bismark'], self.option('sample'))
        if self.option('nucleotide_coverage'):
            cmd += ' --nucleotide_coverage'
        cmd += ' -o {}'.format(self.dir['result'])
        cmd += ' {}'.format(self.option('genome_folder').path)
        cmd += ' -1 {}'.format(self.option('mates1').path)
        cmd += ' -2 {}'.format(self.option('mates2').path)
        if self.option('non_directional'):
            cmd += ' --non_directional'
        runcmd(self, 'run_bismark', cmd)

    def run_deduplicate_bismark(self):
        cmd = '{} --output {} {}'.format(self.program['deduplicate_bismark'], self.dir['result'], self.file['bam'])
        runcmd(self, 'run_deduplicate_bismark', cmd)

    def run_bismark_methylation_extractor(self):
        cmd = '{} -o {}'.format(self.program['bismark_methylation_extractor'], self.dir['result'])
        cmd += ' --parallel {}'.format(self.option('parallel'))
        if self.option('comprehensive'):
            cmd += ' --comprehensive'
        if self.option('bed_graph'):
            cmd += ' --bedGraph'
        if self.option("cx_context"):
            cmd += ' --CX_context'
        if self.option('cytosine_report'):
            cmd += ' --cytosine_report'
            cmd += ' --genome_folder {}'.format(self.option('genome_folder').path)
        cmd += ' {}'.format(self.file['dedup_bam'])
        runcmd(self, 'run_bismark_methylation_extractor', cmd)

    def run_bismark2report(self):
        os.chdir(self.dir['result'])
        cmd = '{} --dir {}'.format(self.program['bismark2report'], self.dir['report'])
        runcmd(self, 'run_bismark2report', cmd)
        os.chdir(self.work_dir)

    def set_output(self):
        pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bismark_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.bismark',
            'instant': False,
            'options': {
                'sample': 'YGXLTHZ1225-4-YG201902876-JJH_L2',
                'genome_folder': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates'
                                 '/Homo_sapiens/Mitochondrion/rCRS',
                'mates1': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/medip/clean_data/YGXLTHZ1225-4'
                          '-YG201902876-JJH_L2.clean.1.fastq',
                'mates2': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/medip/clean_data/YGXLTHZ1225-4'
                          '-YG201902876-JJH_L2.clean.2.fastq',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
