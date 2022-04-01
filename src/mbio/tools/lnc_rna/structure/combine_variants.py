# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class CombineVariantsAgent(Agent):
    '''
    last_modify: 2019.03.12
    '''
    def __init__(self, parent):
        super(CombineVariantsAgent, self).__init__(parent)
        options = [
            {'name': 'vcf_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'genotypemergeoption', 'type': 'string', 'default': 'UNIQUIFY'},
            {'name': 'variant_vcf', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('bcftools')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.bcftools.start()
        self.step.update()

    def step_end(self):
        self.step.bcftools.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} = {}'.format('vcf_list', self.option('vcf_list').path))
        self.logger.debug('{} = {}'.format('ref_fa', self.option('ref_fa').path))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(CombineVariantsAgent, self).end()

class CombineVariantsTool(Tool):
    def __init__(self, config):
        super(CombineVariantsTool, self).__init__(config)
        self.java = 'program/sun_jdk1.8.0/bin/java'
        self.gatk_jar = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/gene-structure/GenomeAnalysisTK.jar')
        self.variant = ','.join([line.strip() for line in open(self.option('vcf_list').path)])
        self.output_vcf = os.path.join(self.work_dir, 'output.vcf')

    def run(self):
        super(CombineVariantsTool, self).run()
        self.run_gatk_combine_variants()
        self.set_output()
        self.end()

    def run_gatk_combine_variants(self):
        cmd = '{} -jar {} -T CombineVariants'.format(self.java, self.gatk_jar)
        cmd += ' -R {}'.format(self.option('ref_fa').path)
        cmd += ' -V {}'.format(self.option('vcf_list').path)
        cmd += ' -o {}'.format(self.output_vcf)
        cmd += ' -genotypeMergeOptions {}'.format(self.option('genotypemergeoption'))
        cmd_name = 'run_gatk_combine_variants'
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
        combine_variants_vcf = os.path.join(self.output_dir, 'combine.variants.vcf')
        if os.path.exists(combine_variants_vcf):
            os.remove(combine_variants_vcf)
        os.link(self.output_vcf, combine_variants_vcf)
        self.logger.info('succeed in linking {} to {}'.format(self.output_vcf, combine_variants_vcf))
        self.option('variant_vcf').set_path(combine_variants_vcf)
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
            'id': 'combine_variants_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.combine_variants',
            'instant': False,
            'options': {
                'vcf_list': '/mnt/ilustre/users/sanger-dev/workspace/20190312/Single_snp_indel_4335_1185/SnpIndel/vcf.list',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/build_idx/output/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fasta',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
