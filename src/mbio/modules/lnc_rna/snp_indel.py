# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan, qindanhua, qinjincheng'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import os
import unittest

class SnpIndelModule(Module):
    '''
    last_modify: 2019.03.13
    '''
    def __init__(self, work_id):
        super(SnpIndelModule, self).__init__(work_id)
        options = [
            {'name': 'method_type', 'type': 'string', 'default': ''},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'bam_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'ref_genome', 'type': 'string', 'default': ''},
            {'name': 'des', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'des_type', 'type': 'string', 'default': ''},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} = {}'.format('method_type', self.option('method_type')))
        self.logger.debug('{} = {}'.format('ref_fa', self.option('ref_fa').path))
        self.logger.debug('{} = {}'.format('bam_list', self.option('bam_list').path))
        self.logger.debug('{} = {}'.format('ref_gtf', self.option('ref_gtf').path))
        self.logger.debug('{} = {}'.format('ref_genome', self.option('ref_genome')))
        self.logger.debug('{} = {}'.format('des', self.option('des').path))
        self.logger.debug('{} = {}'.format('des_type', self.option('des_type')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(SnpIndelModule, self).run()
        if self.option('method_type') == 'samtools':
            self.run_split_bed()
        elif self.option('method_type') == 'gatk':
            self.run_build_idx()

    def run_split_bed(self):
        self.step.add_steps('split_bed')
        self.split_bed = self.add_tool('lnc_rna.structure.split_bed')
        options = {
            'ref_fa': self.option('ref_fa')
        }
        self.split_bed.set_options(options)
        self.split_bed.on('start', self.set_step, {'start': self.step.split_bed})
        self.split_bed.on('end', self.set_step, {'end': self.step.split_bed})
        self.split_bed.on('end', self.run_snptools)
        self.split_bed.run()

    def run_snptools(self):
        for n, bed in enumerate(os.listdir(self.split_bed.output_dir)):
            chr_bed = os.path.join(self.split_bed.output_dir, bed)
            self.step.add_steps('snptools_{}'.format(n))
            snptools = self.add_tool('lnc_rna.structure.snptools')
            snptools.set_options({
                'bam_list': self.option('bam_list'),
                'ref_fa': self.option('ref_fa'),
                'chr_bed': chr_bed,
            })
            snptools.on('start', self.set_step, {'start': getattr(self.step, 'snptools_{}'.format(n))})
            snptools.on('end', self.set_step, {'end': getattr(self.step, 'snptools_{}'.format(n))})
            self.tools.append(snptools)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.run_bcftools)
        else:
            self.on_rely(self.tools, self.run_bcftools)
        for tool in self.tools:
            tool.run()

    def run_bcftools(self):
        vcf_list = os.path.join(self.work_dir, 'vcf.list')
        with open(vcf_list, 'w') as f:
            for tool in self.tools:
                chr_vcf = tool.option('chr_vcf').prop['path']
                if os.path.getsize(chr_vcf) > 0:
                    f.write('{}\n'.format(chr_vcf))
        self.step.add_steps('bcftools')
        self.bcftools = self.add_tool('lnc_rna.structure.bcftools')
        options = {
            'vcf_list': vcf_list
        }
        self.bcftools.set_options(options)
        self.bcftools.on('start', self.set_step, {'start': self.step.bcftools})
        self.bcftools.on('end', self.set_step, {'end': self.step.bcftools})
        self.bcftools.on('end', self.run_annovar)
        self.bcftools.run()

    def run_build_idx(self):
        self.step.add_steps('build_idx')
        self.build_idx = self.add_tool('lnc_rna.structure.build_idx')
        options = {
            'ref_fa': self.option('ref_fa')
        }
        self.build_idx.set_options(options)
        self.build_idx.on('start', self.set_step, {'start': self.step.build_idx})
        self.build_idx.on('end', self.set_step, {'end': self.step.build_idx})
        self.build_idx.on('end', self.run_gatk)
        self.build_idx.run()

    def run_gatk(self):
        for n, line in enumerate(open(self.option('bam_list').path)):
            sample_bam = line.strip()
            self.step.add_steps('gatk_{}'.format(n))
            gatk = self.add_tool('lnc_rna.structure.gatk')
            gatk.set_options({
                'input': sample_bam,
                'ref_fa': self.build_idx.option('reference_fasta')
            })
            gatk.on('start', self.set_step, {'start': getattr(self.step, 'gatk_{}'.format(n))})
            gatk.on('end', self.set_step, {'end': getattr(self.step, 'gatk_{}'.format(n))})
            self.tools.append(gatk)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.run_combine_variants)
        else:
            self.on_rely(self.tools, self.run_combine_variants)
        for tool in self.tools:
            tool.run()

    def run_combine_variants(self):
        vcf_list = os.path.join(self.work_dir, 'vcf.list')
        with open(vcf_list, 'w') as f:
            for tool in self.tools:
                sample_vcf = tool.option('sample_vcf').prop['path']
                if os.path.getsize(sample_vcf) > 0:
                    f.write('{}\n'.format(sample_vcf))
        self.step.add_steps('combine_variants')
        self.combine_variants = self.add_tool('lnc_rna.structure.combine_variants')
        options = {
            'vcf_list': vcf_list,
            'ref_fa': self.build_idx.option('reference_fasta')
        }
        self.combine_variants.set_options(options)
        self.combine_variants.on('start', self.set_step, {'start': self.step.combine_variants})
        self.combine_variants.on('end', self.set_step, {'end': self.step.combine_variants})
        self.combine_variants.on('end', self.run_annovar)
        self.combine_variants.run()

    def run_annovar(self):
        self.step.add_steps('annovar')
        self.annovar = self.add_tool('lnc_rna.structure.annovar')
        options = {
            'ref_gtf': self.option('ref_gtf'),
            'ref_fa': self.option('ref_fa'),
            'ref_genome': self.option('ref_genome'),
            'des': self.option('des'),
            'des_type': self.option('des_type'),
        }
        if self.option('method_type') == 'samtools':
            options.update({'variant_vcf': self.bcftools.option('variant_vcf')})
        elif self.option('method_type') == 'gatk':
            options.update({'variant_vcf': self.combine_variants.option('variant_vcf')})
        self.annovar.set_options(options)
        self.annovar.on('start', self.set_step, {'start': self.step.annovar})
        self.annovar.on('end', self.set_step, {'end': self.step.annovar})
        self.annovar.on('end', self.run_snp_stat)
        self.annovar.run()

    def run_snp_stat(self):
        self.step.add_steps('snp_stat')
        self.snp_stat = self.add_tool('lnc_rna.structure.snp_stat')
        options = {
            'snp_detail': self.annovar.option('snp_detail'),
            'snp_annotation': self.annovar.option('snp_annotation'),
        }
        self.snp_stat.set_options(options)
        self.snp_stat.on('start', self.set_step, {'start': self.step.snp_stat})
        self.snp_stat.on('end', self.set_step, {'end': self.step.snp_stat})
        self.snp_stat.on('end', self.set_output)
        self.snp_stat.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for basename in os.listdir(self.snp_stat.output_dir):
            source = os.path.join(self.snp_stat.output_dir, basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(SnpIndelModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    # def test_samtools(self):
    #     import random
    #     from mbio.workflows.single import SingleWorkflow
    #     from biocluster.wsheet import Sheet
    #     data = {
    #         'id': 'snp_indel_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
    #         'type': 'module',
    #         'name': 'lnc_rna.snp_indel',
    #         'instant': False,
    #         'options': {
    #             'method_type': 'samtools',
    #             'bam_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/test_data/bam_list.txt',
    #             'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
    #             'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
    #             'ref_genome': 'Homo_sapiens',
    #             'des': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
    #             'des_type': 'type1',
    #         }
    #     }
    #     wsheet = Sheet(data=data)
    #     wf = SingleWorkflow(wsheet)
    #     wf.run()

    def test_gatk(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'snp_indel_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.snp_indel',
            'instant': False,
            'options': {
                'method_type': 'gatk',
                'bam_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/test_data/bam_list.txt',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_genome': 'Homo_sapiens',
                'des': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                'des_type': 'type1',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()