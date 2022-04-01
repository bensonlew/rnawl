# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class AnnovarAgent(Agent):
    '''
    last_modify: 2019.03.06
    '''
    def __init__(self, parent):
        super(AnnovarAgent, self).__init__(parent)
        options = [
            {'name': 'variant_vcf', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'ref_genome', 'type': 'string', 'default': ''},
            {'name': 'prefix', 'type': 'string', 'default': 'output'},
            {'name': 'dbtype', 'type': 'string', 'default': 'refGene'},
            {'name': 'des', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'des_type', 'type': 'string', 'default': ''},
            {'name': 'snp_detail', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'snp_annotation', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('annovar')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.annovar.start()
        self.step.update()

    def step_end(self):
        self.step.annovar.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('variant_vcf').is_set:
            self.logger.debug('{} = {}'.format('variant_vcf', self.option('variant_vcf').path))
            self.infile_size = os.path.getsize(self.option('variant_vcf').path)
            self.fasta_size = os.path.getsize(self.option('ref_fa').path)
        else:
            raise OptionError('VCF file must be provided')
        if self.option('ref_gtf').is_set:
            self.logger.debug('{} = {}'.format('ref_gtf', self.option('ref_gtf').path))
        if self.option('ref_fa').is_set:
            self.logger.debug('{} = {}'.format('ref_fa', self.option('ref_fa').path))
        self.logger.debug('{} = {}'.format('ref_genome', self.option('ref_genome')))
        self.logger.debug('{} = {}'.format('prefix', self.option('prefix')))
        self.logger.debug('{} = {}'.format('dbtype', self.option('dbtype')))
        if self.option('des').is_set:
            self.logger.debug('{} = {}'.format('des', self.option('des').path))
        self.logger.debug('{} = {}'.format('des_type', self.option('des_type')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        mem = max(int(float(self.fasta_size)/1024 ** 3 + 10) , int(float(self.infile_size) / 1024 ** 3 * 8 + 16))
        self._memory = '{}G'.format(mem)

    def end(self):
        super(AnnovarAgent, self).end()

class AnnovarTool(Tool):
    def __init__(self, config):
        super(AnnovarTool, self).__init__(config)
        self.gtftogenepred = 'bioinfo/gene-structure/annovar/gtfToGenePred'
        self.python = 'program/Python/bin/python'
        self.perl = 'program/perl-5.24.0/bin/perl'
        self.clean_vcf_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/clean_vcf.py')
        self.genepred2refgene_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/genepred2refgene.py')
        self.convert2annovar_pl = os.path.join(
            self.config.SOFTWARE_DIR, 'bioinfo/gene-structure/annovar/convert2annovar.pl'
        )
        self.retrieve_seq_from_fasta_pl = os.path.join(
            self.config.SOFTWARE_DIR, 'bioinfo/gene-structure/annovar/retrieve_seq_from_fasta.pl'
        )
        self.annotate_variation_pl = os.path.join(
            self.config.SOFTWARE_DIR, 'bioinfo/gene-structure/annovar/annotate_variation.pl'
        )
        self.snp_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/snp_annotation.py')
        self.get_id_name_des_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/get_id_name_des.py')
        self.add_name_des_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/add_name_des.py')

        self.vcf = os.path.join(self.work_dir, 'input.vcf')
        self.genepred = os.path.join(self.work_dir, '{}_genePred.txt'.format(self.option('ref_genome')))
        self.database = os.path.join(self.work_dir, 'database')
        if not os.path.isdir(self.database):
            os.mkdir(self.database)
        self.refgene = os.path.join(self.database, '{}_refGene.txt'.format(self.option('ref_genome')))
        self.input_fa = os.path.join(self.database, '{}_refGeneMrna.fa'.format(self.option('ref_genome')))
        self.avinput = os.path.join(self.work_dir, '{}.avinput'.format(self.option('prefix')))
        self.exonic_variant_function = os.path.join(
            self.work_dir, '{}.exonic_variant_function'.format(self.option('prefix'))
        )
        self.variant_function = os.path.join(self.work_dir, '{}.variant_function'.format(self.option('prefix')))
        self.detail_xls = os.path.join(self.work_dir, 'snp_detail.xls')
        self.id_name_des_txt = os.path.join(self.work_dir, 'id_name_des.txt')
        self.annotation_xls = os.path.join(self.work_dir, 'snp_annotation.xls')

    def run(self):
        super(AnnovarTool, self).run()
        self.run_clean_vcf()
        self.run_convert2annovar()
        self.run_gtftogenepred()
        self.run_genepred2refgene()
        self.run_retrieve_seq_from_fasta()
        self.run_annotate_variation()
        self.run_snp_annotation()
        self.run_get_id_name_des()
        self.run_add_name_des()
        self.set_output()
        self.end()

    def run_clean_vcf(self):
        cmd = '{} {}'.format(self.python, self.clean_vcf_py)
        cmd += ' -i {}'.format(self.option('variant_vcf').path)
        cmd += ' -o {}'.format(self.vcf)
        cmd_name = 'run_clean_vcf'
        self.run_code(cmd_name, cmd)

    def run_convert2annovar(self):
        cmd = '{} {}'.format(self.perl, self.convert2annovar_pl)
        cmd += ' --format {}'.format('vcf4old')
        cmd += ' --allallele'
        cmd += ' --outfile {}'.format(self.avinput)
        cmd += ' {}'.format(self.vcf)
        cmd_name = 'run_convert2annovar'
        self.run_code(cmd_name, cmd)

    def run_gtftogenepred(self):
        cmd = '{} -genePredExt -allErrors {} {}'.format(
            self.gtftogenepred,
            self.option('ref_gtf').path,
            self.genepred
        )
        cmd_name = 'run_gtftogenepred'
        self.run_code(cmd_name, cmd)

    def run_genepred2refgene(self):
        cmd = '{} {}'.format(self.python, self.genepred2refgene_py)
        cmd += ' -i {}'.format(self.genepred)
        cmd += ' -o {}'.format(self.refgene)
        cmd_name = 'run_genepred2refgene'
        self.run_code(cmd_name, cmd)

    def run_retrieve_seq_from_fasta(self):
        cmd = '{} {} {}'.format(self.perl, self.retrieve_seq_from_fasta_pl, self.refgene)
        cmd += ' --format {}'.format('refGene')
        cmd += ' --outfile {}'.format(self.input_fa)
        cmd += ' --seqfile {}'.format(self.option('ref_fa').path)
        cmd_name = 'run_retrieve_seq_from_fasta'
        self.run_code(cmd_name, cmd)

    def run_annotate_variation(self):
        cmd = '{} {} {} {}'.format(self.perl, self.annotate_variation_pl, self.avinput, self.database)
        cmd += ' --outfile {}'.format(self.option('prefix'))
        cmd += ' --dbtype {}'.format(self.option('dbtype'))
        cmd += ' --buildver {}'.format(self.option('ref_genome'))
        cmd_name = 'run_annotate_variation'
        self.run_code(cmd_name, cmd)

    def run_snp_annotation(self):
        cmd = '{} {}'.format(self.python, self.snp_annotation_py)
        cmd += ' -i {}'.format(self.vcf)
        cmd += ' -e {}'.format(self.exonic_variant_function)
        cmd += ' -v {}'.format(self.variant_function)
        cmd += ' -o {}'.format(self.detail_xls)
        cmd_name = 'run_snp_annotation'
        self.run_code(cmd_name, cmd)

    def run_get_id_name_des(self):
        cmd = '{} {}'.format(self.python, self.get_id_name_des_py)
        cmd += ' -i {}'.format(self.option('des').path)
        cmd += ' -t {}'.format(self.option('des_type'))
        cmd += ' -o {}'.format(self.id_name_des_txt)
        cmd_name = 'run_get_id_name_des'
        self.run_code(cmd_name, cmd)

    def run_add_name_des(self):
        cmd = '{} {}'.format(self.python, self.add_name_des_py)
        cmd += ' -i {}'.format(self.detail_xls)
        cmd += ' -a {}'.format(self.id_name_des_txt)
        cmd += ' -o {}'.format(self.annotation_xls)
        cmd_name = 'run_add_name_des'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        snp_detail = os.path.join(self.output_dir, 'snp_detail.xls')
        if os.path.exists(snp_detail):
            os.remove(snp_detail)
        os.link(self.detail_xls, snp_detail)
        self.logger.info('succeed in linking {} to {}'.format(self.detail_xls, snp_detail))
        self.option('snp_detail').set_path(snp_detail)
        snp_annotation = os.path.join(self.output_dir, 'snp_annotation.xls')
        if os.path.exists(snp_annotation):
            os.remove(snp_annotation)
        os.link(self.annotation_xls, snp_annotation)
        self.logger.info('succeed in linking {} to {}'.format(self.annotation_xls, snp_annotation))
        self.option('snp_annotation').set_path(snp_annotation)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_hsa(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annovar_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.annovar',
            'instant': False,
            'options': {
                'variant_vcf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/bcftools/output/pop.noid.vcf',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
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
