# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class CircfinderAgent(Agent):
    '''
    last_modify: 2019.8.22
    '''

    def __init__(self, parent):
        super(CircfinderAgent, self).__init__(parent)
        options = [
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'fq1', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'fq2', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'star_index', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 's_filteredJunctions', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'circ_finder', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'circ_finder_filter', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'star_sam', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'junction_reads', 'type': 'int', 'default': 2},
            {'name': 'circrna_length', 'type': 'int', 'default': 100000}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '80G'

    def end(self):
        super(CircfinderAgent, self).end()


class CircfinderTool(Tool):
    def __init__(self, config):
        super(CircfinderTool, self).__init__(config)
        self.program = {
            'star': 'bioinfo/rna/star-2.5/bin/Linux_x86_64/',
            'python': 'miniconda2/bin/python',
            'perl': 'program/perl-5.24.0/bin/perl'

        }
        self.file = {
            'star_index': self.output_dir,
            'circfinder_dir': self.work_dir,
            's_filteredJunctions': os.path.join(self.work_dir, 's_filteredJunctions.bed'),
            'circ_finder': os.path.join(self.output_dir, '{}_circ_finder.bed'.format(self.option('sample'))),
            'circ_finder_filter': os.path.join(self.output_dir,
                                               '{}_circ_finder_filter.txt'.format(self.option('sample'))),
            'star_sam': os.path.join(self.work_dir, 'Aligned.out.sam')

        }
        self.script = {
            'circ_finder': os.path.join(self.config.SOFTWARE_DIR,
                                        'bioinfo/circ_rna/circRNA_finder-master/postProcessStarAlignment.pl')
        }

    def run(self):
        super(CircfinderTool, self).run()
        self.star_index()
        self.star_aln_pe()
        self.star_aln_pe_bam()
        self.run_circ_finder()
        self.run_rename()
        self.filter()
        self.set_output()
        self.end()

    def star_index(self):
        """
        step1:第一步建索引；用star建立参考基因组的索引，当用户不上传参考基因组时，该步骤省略，直接调用已有的序列文件
        genomeDir为用于存放第一步建立的参考基因组索引的路径

        """
        cmd = "{}STAR --runMode genomeGenerate --limitGenomeGenerateRAM 50000000000 --genomeDir {} --genomeFastaFiles {} --sjdbGTFfile {} --runThreadN 10".format(
            self.program['star'], self.file['star_index'], self.option('genome').path, self.option('annotate').path)
        cmd_name = 'star_index'
        runcmd(self, cmd_name, cmd)

    def star_aln_pe(self):
        """
        step2:第二步比对；用star进行双端序列的比对
        """
        if self.option('fq2').is_set:
            cmd = "{}STAR --runThreadN 10 --genomeDir {} --readFilesIn {} {}".format(self.program['star'],
                                                                                     self.file['star_index'],
                                                                                     self.option("fq1").prop["path"],
                                                                                     self.option("fq2").prop["path"])
            cmd_name = 'star_aln1_pe'
            runcmd(self, cmd_name, cmd)
        else:
            cmd = "{}STAR --runThreadN 10 --genomeDir {} --readFilesIn {}".format(self.program['star'],
                                                                                  self.file['star_index'],
                                                                                  self.option("fq1").prop["path"])
            cmd_name = 'star_aln1_pe'
            runcmd(self, cmd_name, cmd)

    def star_aln_pe_bam(self):
        """
        step4：第四步，最终比对
        """
        if self.option('fq2').is_set:
            cmd = "{}STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --genomeDir {} --readFilesIn {} {}".format(
                self.program['star'], self.file['star_index'], self.option("fq1").prop["path"],
                self.option("fq2").prop["path"])
        else:
            cmd = "{}STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --genomeDir {} --readFilesIn {}".format(
                self.program['star'], self.file['star_index'], self.option("fq1").prop["path"])
        cmd += ' --chimSegmentMin 10'
        cmd += ' --chimScoreMin 1'
        cmd += ' --alignIntronMax 500000'
        cmd += ' --outFilterMismatchNmax 4'
        cmd += ' --alignTranscriptsPerReadNmax 100000'
        cmd += ' --chimOutType Junctions SeparateSAMolFd --outFilterMultimapNmax 2'
        cmd_name = 'star_aln2_pe_bam'
        runcmd(self, cmd_name, cmd)

    def run_circ_finder(self):
        cmd = '{} {}'.format(self.program['perl'], self.script['circ_finder'])
        cmd_name = 'run_circ_finder'
        runcmd(self, cmd_name, cmd)

    def run_rename(self):
        os.rename(self.file['s_filteredJunctions'],
                  self.file['circ_finder'])
        # cmd = 'mv {} {}'.format(self.file['s_filteredJunctions'], self.file['circ_finder'])
        # cmd_name = 'run_rename'
        # runcmd(self, cmd_name,cmd,shell = True)

    def filter(self):
        f1 = pd.read_table(self.file['circ_finder'], header=None)
        f2 = f1[[0, 1, 2, 4, 5]]
        f3 = f2.rename(columns={0: 'chr', 1: 'circRNA_start', 2: 'circRNA_end', 4: 'junction_reads', 5: 'strand'})
        f4 = f3[f3['circRNA_end'] - f3['circRNA_start'] < self.option('circrna_length')]
        f5 = f4[f4['junction_reads'] >= self.option('junction_reads')]
        f5['circRNA_start'] = f5['circRNA_start'] + 1
        index = f5.index.tolist()
        f5['circRNA_name'] = ['{}_{}_{}'.format(f5['chr'][i], f5['circRNA_start'][i], f5['circRNA_end'][i]) for i in
                              index]
        f5.sort_values(by=['chr', 'circRNA_start', 'circRNA_end'], ascending=True, inplace=True)
        f5.to_csv(self.file['circ_finder_filter'], index=False, sep='\t')

    def set_output(self):
        self.option('circ_finder_filter').set_path(self.file['circ_finder_filter'])
        self.option('star_sam').set_path(self.file['star_sam'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'circ_finder_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.circfinder',
            'instant': False,
            'options': {

                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.dna.toplevel.fa',
                # 'annotate_txt': '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/circRNA/circ_predict_software/test_data/test_pre/test_anno/hg19_ref_all.txt',
                'annotate': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.96.gtf',
                'fq1': '/mnt/ilustre/users/isanger/workspace/20190213/Single_LncrnaQc_7789/LncrnaQc/output/sickle_dir/Con1_sickle_l.fastq',
                'fq2': '/mnt/ilustre/users/isanger/workspace/20190213/Single_LncrnaQc_7789/LncrnaQc/output/sickle_dir/Con1_sickle_r.fastq',
                'sample': 'zjx'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
