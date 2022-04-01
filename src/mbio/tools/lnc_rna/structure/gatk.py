# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class GatkAgent(Agent):
    '''
    last_modify: 2019.05.10
    '''
    def __init__(self, parent):
        super(GatkAgent, self).__init__(parent)
        options = [
            {'name': 'input', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'validation_stringency', 'type': 'string', 'default': 'SILENT'},
            {'name': 'sort_order', 'type': 'string', 'default': 'coordinate'},
            {'name': 'rglb', 'type': 'string', 'default': 'lib1'},
            {'name': 'rgpl', 'type': 'string', 'default': 'illumina'},
            {'name': 'rgpu', 'type': 'string', 'default': 'unit1'},
            {'name': 'create_index', 'type': 'bool', 'default': True},
            {'name': 'ref_fa', 'type': 'infile', 'format':'lnc_rna.fasta'},
            {'name': 'unsafe', 'type': 'string', 'default': 'ALLOW_N_CIGAR_READS'},
            {'name': 'read_filter', 'type': 'string', 'default': 'ReassignOneMappingQuality'},
            {'name': 'cpu', 'type': 'int', 'default': 8},
            {'name': 'cluster', 'type': 'int', 'default': 3},
            {'name': 'window', 'type': 'int', 'default': 35},
            {'name': 'max_FS', 'type': 'float', 'default': 30.0},
            {'name': 'min_QD', 'type': 'float', 'default': 2.0},
            {'name': 'sample_vcf', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('gatk')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.gatk.start()
        self.step.update()

    def step_end(self):
        self.step.gatk.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('input').is_set:
            self.logger.debug('{} = {}'.format('input', self.option('input').prop['path']))
        self.logger.debug('{} = {}'.format('validation_stringency', self.option('validation_stringency')))
        self.logger.debug('{} = {}'.format('sort_order', self.option('sort_order')))
        self.logger.debug('{} = {}'.format('rglb', self.option('rglb')))
        self.logger.debug('{} = {}'.format('rgpl', self.option('rgpl')))
        self.logger.debug('{} = {}'.format('rgpu', self.option('rgpu')))
        self.logger.debug('{} = {}'.format('create_index', self.option('create_index')))
        if self.option('ref_fa').is_set:
            self.logger.debug('{} = {}'.format('ref_fa', self.option('ref_fa').prop['path']))
            self.infile_size = os.path.getsize(self.option('ref_fa').prop['path'])
        self.logger.debug('{} = {}'.format('unsafe', self.option('unsafe')))
        self.logger.debug('{} = {}'.format('read_filter', self.option('read_filter')))
        self.logger.debug('{} = {}'.format('cpu', self.option('cpu')))
        self.logger.debug('{} = {}'.format('cluster', self.option('cluster')))
        self.logger.debug('{} = {}'.format('window', self.option('window')))
        self.logger.debug('{} = {}'.format('max_FS', self.option('max_FS')))
        self.logger.debug('{} = {}'.format('min_QD', self.option('min_QD')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('cpu')
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 8 + 24))

    def end(self):
        super(GatkAgent, self).end()

class GatkTool(Tool):
    def __init__(self, config):
        super(GatkTool, self).__init__(config)
        self.java = 'program/sun_jdk1.8.0/bin/java'
        self.picard_jar = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/gene-structure/picard.jar')
        self.gatk_jar = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/gene-structure/GenomeAnalysisTK.jar')
        self.sample = os.path.basename(self.option('input').path[:-4])
        self.add_sorted_bam = os.path.join(self.work_dir, 'add_sorted.bam')
        self.marked_duplicates_bam = os.path.join(self.work_dir, 'marked_duplicates.bam')
        self.marked_dup_metrics_txt = os.path.join(self.work_dir, 'marked_dup_metrics.txt')
        self.output_bam = os.path.join(self.work_dir, 'output.bam')
        self.output_vcf = os.path.join(self.work_dir, 'output.vcf')
        self.filter_vcf = os.path.join(self.work_dir, 'filter.vcf')

    def run(self):
        super(GatkTool, self).run()
        self.run_picard_add_or_replace_read_groups()
        self.run_picard_mark_duplicates()
        self.run_picard_build_bam_index()
        self.run_gatk_split_n_cigar_reads()
        self.run_gatk_haplotype_caller()
        self.run_gatk_variant_filtration()
        self.set_output()
        self.end()

    def run_picard_add_or_replace_read_groups(self):
        cmd = '{} -jar {} AddOrReplaceReadGroups'.format(self.java, self.picard_jar)
        cmd += ' VALIDATION_STRINGENCY={}'.format(self.option('validation_stringency'))
        cmd += ' I={}'.format(self.option('input').path)
        cmd += ' O={}'.format(self.add_sorted_bam)
        cmd += ' SORT_ORDER={}'.format(self.option('sort_order'))
        cmd += ' RGLB={}'.format(self.option('rglb'))
        cmd += ' RGPL={}'.format(self.option('rgpl'))
        cmd += ' RGPU={}'.format(self.option('rgpu'))
        cmd += ' RGSM={}'.format(self.sample)
        cmd_name = 'run_picard_add_or_replace_read_groups'
        self.run_code(cmd_name, cmd)

    def run_picard_mark_duplicates(self):
        cmd = '{} -jar {} MarkDuplicates'.format(self.java, self.picard_jar)
        cmd += ' VALIDATION_STRINGENCY={}'.format(self.option('validation_stringency'))
        cmd += ' CREATE_INDEX={}'.format(self.option('create_index'))
        cmd += ' I={}'.format(self.add_sorted_bam)
        cmd += ' O={}'.format(self.marked_duplicates_bam)
        cmd += ' M={}'.format(self.marked_dup_metrics_txt)
        cmd_name = 'run_picard_mark_duplicates'
        self.run_code(cmd_name, cmd)

    def run_picard_build_bam_index(self):
        cmd = '{} -jar {} BuildBamIndex'.format(self.java, self.picard_jar)
        cmd += ' I={}'.format(self.marked_duplicates_bam)
        cmd_name = 'run_picard_build_bam_index'
        self.run_code(cmd_name, cmd)

    def run_gatk_split_n_cigar_reads(self):
        cmd = '{} -jar {} -T SplitNCigarReads'.format(self.java, self.gatk_jar)
        cmd += ' -R {}'.format(self.option('ref_fa').path)
        cmd += ' -I {}'.format(self.marked_duplicates_bam)
        cmd += ' -o {}'.format(self.output_bam)
        cmd += ' -U {}'.format(self.option('unsafe'))
        cmd += ' -rf {}'.format(self.option('read_filter'))
        cmd_name = 'run_gatk_split_n_cigar_reads'
        self.run_code(cmd_name, cmd)

    def run_gatk_haplotype_caller(self):
        cmd = '{} -jar {} -T HaplotypeCaller'.format(self.java, self.gatk_jar)
        cmd += ' -R {}'.format(self.option('ref_fa').path)
        cmd += ' -I {}'.format(self.output_bam)
        cmd += ' -o {}'.format(self.output_vcf)
        cmd_name = 'run_gatk_haplotype_caller'
        self.run_code(cmd_name, cmd)

    def run_gatk_variant_filtration(self):
        cmd = '{} -jar {} -T VariantFiltration'.format(self.java, self.gatk_jar)
        cmd += ' -R {}'.format(self.option('ref_fa').path)
        cmd += ' -V {}'.format(self.output_vcf)
        cmd += ' -o {}'.format(self.filter_vcf)
        cmd += ' -cluster {}'.format(self.option('cluster'))
        cmd += ' -window {}'.format(self.option('window'))
        cmd += ' -filter "FS > {}" -filterName FS'.format(self.option('max_FS'))
        cmd += ' -filter "QD < {}" -filterName QD'.format(self.option('min_QD'))
        cmd_name = 'run_gatk_variant_filtration'
        self.run_code(cmd_name, cmd, shell=True)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
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
        filter_vcf = os.path.join(self.output_dir, os.path.basename(self.filter_vcf))
        if os.path.exists(filter_vcf):
            os.remove(filter_vcf)
        os.link(self.filter_vcf, filter_vcf)
        self.logger.info('succeed in linking {} to {}'.format(self.filter_vcf, filter_vcf))
        self.option('sample_vcf').set_path(filter_vcf)
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
            'id': 'gatk_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.gatk',
            'instant': False,
            'options': {
                'input': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/gatk/Con1.bam',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/build_idx/output/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fasta',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()