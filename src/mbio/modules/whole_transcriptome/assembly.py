# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.module import Module


class AssemblyModule(Module):
    def __init__(self, work_id):
        super(AssemblyModule, self).__init__(work_id)
        options = [
            {'name': 'bamlist', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'strand_specific', 'type': 'bool', 'default': True},
            {'name': 'strand_dir', 'type': 'string', 'default': 'RF'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(AssemblyModule, self).run()
        self.run_stringtie()

    def run_stringtie(self):
        for line in open(self.option('bamlist').path):
            stringtie = self.add_tool('whole_transcriptome.assembly.stringtie')
            input_bam = line.strip()
            opts = {
                'input_bam': input_bam,
                'strand_specific': self.option('strand_specific'),
                'strand_dir': self.option('strand_dir'),
                'guide_gff': self.option('ref_gtf')
            }
            stringtie.set_options(opts)
            self.tools.append(stringtie)
        else:
            self.on_rely(self.tools, self.run_stringtie_merge)
        for tool in self.tools:
            tool.run()

    def run_stringtie_merge(self):
        self.stringtie_merge = self.add_tool('whole_transcriptome.assembly.stringtie_merge')
        gtf_list = os.path.join(self.work_dir, 'gtf.list')
        open(gtf_list, 'w').writelines('{}\n'.format(tool.option('out_gtf').path) for tool in self.tools)
        opts = {
            'gtf_list': gtf_list,
            'guide_gff': self.option('ref_gtf'),
        }
        self.stringtie_merge.set_options(opts)
        self.stringtie_merge.on('end', self.run_trans_build)
        self.stringtie_merge.run()

    def run_trans_build(self):
        self.trans_build = self.add_tool('whole_transcriptome.assembly.trans_build')
        opts = {
            'ref_gtf': self.option('ref_gtf'),
            'seq_path': self.option('ref_fa'),
            'merged_gtf': self.stringtie_merge.option('out_gtf'),
        }
        self.trans_build.set_options(opts)
        self.trans_build.on('end', self.set_output)
        self.trans_build.run()

    def set_output(self):
        for file_name in os.listdir(self.trans_build.output_dir):
            source = os.path.join(self.trans_build.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        else:
            self.end()

    def end(self):
        super(AssemblyModule, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'assembly_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.assembly',
            'instant': False,
            'options': {
                'bamlist': '/mnt/lustre/users/sanger/workspace/20191230/Longrna_sanger_231556/RnaseqMapping/output/bamlist',
                'strand_specific': True,
                'strand_dir': 'RF',
                'ref_gtf': '/mnt/lustre/users/sanger/workspace/20191230/Longrna_sanger_231556/FileCheck/Mus_musculus.GRCm38.96.gtf',
                'ref_fa': '/mnt/lustre/users/sanger/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/dna/Mus_musculus.GRCm38.dna.toplevel.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
