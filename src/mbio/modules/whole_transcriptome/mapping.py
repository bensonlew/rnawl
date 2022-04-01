# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.module import Module
from mbio.packages.whole_transcriptome.utils import read_fastq_dir

PROGRAM = ('bowtie2', 'hisat2')


class MappingModule(Module):
    def __init__(self, work_id):
        super(MappingModule, self).__init__(work_id)
        options = [
            {'name': 'program', 'type': 'string', 'default': PROGRAM[0]},
            {'name': 'index', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'strand_specific', 'type': 'bool', 'default': False},
            {'name': 'strand_dir', 'type': 'string', 'default': 'RF'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} => {}'.format(k, v.value))

    def run(self):
        super(MappingModule, self).run()
        if self.option('program') == 'bowtie2':
            self.run_bowtie2()
        elif self.option('program') == 'hisat2':
            self.run_hisat2()

    def run_bowtie2(self):
        is_se, fastqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        bt2_idx = self.option('index').path
        for sample_name, fastqs in fastqs_dict.items():
            bowtie2 = self.add_tool('whole_transcriptome.bowtie2')
            opts = {'bt2_idx': bt2_idx, 'is_se': is_se, 'sample_name': sample_name}
            if is_se:
                opts['r'] = fastqs[0]
            else:
                opts['m1'] = fastqs[0]
                opts['m2'] = fastqs[1]
            bowtie2.set_options(opts)
            self.tools.append(bowtie2)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def run_hisat2(self):
        is_se, fastqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        ht2_idx = self.option('index').path
        for sample_name, fastqs in fastqs_dict.items():
            hisat2 = self.add_tool('whole_transcriptome.hisat2')
            opts = {'ht2_idx': ht2_idx, 'is_se': is_se, 'sample_name': sample_name}
            if is_se:
                opts['r'] = fastqs[0]
            else:
                opts['m1'] = fastqs[0]
                opts['m2'] = fastqs[1]
            hisat2.set_options(opts)
            self.tools.append(hisat2)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        bam_dir = os.path.join(self.output_dir, 'map_data')
        if os.path.isdir(bam_dir):
            shutil.rmtree(bam_dir)
        os.mkdir(bam_dir)
        for tool in self.tools:
            source = tool.option('bam').path
            link_name = os.path.join(bam_dir, os.path.basename(source))
            os.link(source, link_name)
        self.end()

    def end(self):
        super(MappingModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the module. Just run this script to do test.
    """

    def test_bowtie2(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mapping_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.mapping',
            'instant': False,
            'options': {
                'program': 'bowtie2',
                'index': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens'
                         '/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/medip/clean_data',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_hisat2(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mapping_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.mapping',
            'instant': False,
            'options': {
                'program': 'hisat2',
                'index': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/bacilli'
                         '/Geobacillus_thermodenitrificans/GCF_000015745.1_ASM1574v1/dna/GCF_000015745'
                         '.1_ASM1574v1_genomic',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/data/Geobacillus_thermodenitrificans'
                             '/clean_data',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_hisat2_o2(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mapping_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.mapping',
            'instant': False,
            'options': {
                'program': 'hisat2',
                'index': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/bacilli'
                         '/Geobacillus_thermodenitrificans/GCF_000015745.1_ASM1574v1/dna/GCF_000015745'
                         '.1_ASM1574v1_genomic',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/data'
                             '/Geobacillus_thermodenitrificans_O2/clean_data',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_hisat2_o2')])
    unittest.TextTestRunner(verbosity=2).run(suite)
