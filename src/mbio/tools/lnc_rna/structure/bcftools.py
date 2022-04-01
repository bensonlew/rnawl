# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class BcftoolsAgent(Agent):
    '''
    last_modify: 2019.03.01
    '''
    def __init__(self, parent):
        super(BcftoolsAgent, self).__init__(parent)
        options = [
            {'name': 'vcf_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'output_type', 'type': 'string', 'default': 'v'},
            {'name': 'max_mem', 'type': 'string', 'default': '20G'},
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
        if self.option('vcf_list').is_set:
            self.logger.debug('{} = {}'.format('vcf_list', self.option('vcf_list').prop['path']))
        else:
            raise OptionError('VCF list file must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = self.option('max_mem')

    def end(self):
        super(BcftoolsAgent, self).end()

class BcftoolsTool(Tool):
    def __init__(self, config):
        super(BcftoolsTool, self).__init__(config)
        self.bcftools = 'bioinfo/seq/bcftools-1.7/bcftools'
        self.pop_nosort_vcf = os.path.join(self.work_dir, 'pop.nosort.vcf')
        self.pop_noid_vcf = os.path.join(self.work_dir, 'pop.noid.vcf')

    def run(self):
        super(BcftoolsTool, self).run()
        self.run_bcftools_concat()
        self.run_bcftools_sort()
        self.set_output()
        self.end()

    def run_bcftools_concat(self):
        cmd = '{} concat'.format(self.bcftools)
        cmd += ' -f {}'.format(self.option('vcf_list').prop['path'])
        cmd += ' -o {}'.format(self.pop_nosort_vcf)
        cmd += ' -O {}'.format(self.option('output_type'))
        cmd_name = 'run_bcftools_concat'
        self.run_code(cmd_name, cmd)

    def run_bcftools_sort(self):
        cmd = '{} sort {}'.format(self.bcftools, self.pop_nosort_vcf)
        cmd += ' -m {}'.format(self.option('max_mem'))
        cmd += ' -o {}'.format(self.pop_noid_vcf)
        cmd += ' -O {}'.format(self.option('output_type'))
        cmd_name = 'run_bcftools_sort'
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
        noid_vcf = os.path.join(self.output_dir, 'pop.variant.vcf')
        if os.path.exists(noid_vcf):
            os.remove(noid_vcf)
        os.link(self.pop_noid_vcf, noid_vcf)
        self.logger.info('succeed in linking {} to {}'.format(self.pop_noid_vcf, noid_vcf))
        self.option('variant_vcf').set_path(noid_vcf)
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
            'id': 'bcftools_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.bcftools',
            'instant': False,
            'options': {
                'vcf_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/snptools/output/vcf.list',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
