# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.packages.ref_rna_v3.functions import toolfuncdeco, runcmd


class RpkmSaturationAgent(Agent):
    '''
    last_modify: 2019.07.05
    '''

    def __init__(self, parent):
        super(RpkmSaturationAgent, self).__init__(parent)
        options = [
            {'name': 'bed', 'type': 'infile', 'format': 'ref_rna_v2.bed'},  # bed格式文件
            {'name': 'bam', 'type': 'infile', 'format': 'align.bwa.bam'},  # bam格式文件，排序过的
            {'name': 'quality', 'type': 'int', 'default': 30},  # Minimum mapping quality
            {'name': 'low_bound', 'type': 'int', 'default': 5},  # Sampling starts from this percentile
            {'name': 'up_bound', 'type': 'int', 'default': 100},  # Sampling ends at this percentile
            {'name': 'step', 'type': 'int', 'default': 5},  # Sampling frequency
            {'name': 'rpkm_cutof', 'type': 'float', 'default': 0.01}
            # Transcripts with RPKM smaller than this number will be ignored
        ]
        self.add_option(options)
        self._memory_increase_step = 100
        self.step.add_steps('satur')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.satur.start()
        self.step.update()

    def step_end(self):
        self.step.satur.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        if not self.option('bam').is_set:
            raise OptionError('请传入比对结果bam格式文件', code='35601103')
        if not self.option('bed').is_set:
            raise OptionError('请传入bed格式文件', code='35601104')
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '40G'
        bam_size = float(os.path.getsize(self.option('bam').prop['path'])) / 1024 ** 3
        if bam_size > 5:
            self._memory = '{}G'.format(int((bam_size - 5)*5 + 40))

    @toolfuncdeco
    def end(self):
        super(RpkmSaturationAgent, self).end()


class RpkmSaturationTool(Tool):
    def __init__(self, config):
        super(RpkmSaturationTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin'))
        self.prefix = '{}_{}'.format(
            os.path.join(self.work_dir, 'satur'), os.path.basename(self.option('bam').path)[:-4]
        )
        self.program = {
            'python': 'miniconda2/bin/python',
            'perl': 'program/perl/perls/perl-5.24.0/bin/perl'
        }
        self.script = {
            'RPKM_saturation': os.path.join(self.config.SOFTWARE_DIR,
                                            'bioinfo/ref_rna_v2/miniconda2/bin/RPKM_saturation.py'),
            'saturation2plot': os.path.join(self.config.SOFTWARE_DIR,
                                            'bioinfo/plot/scripts/saturation2plot.pl')
        }
        self.file = {
            'xls': '{}.eRPKM.xls'.format(self.prefix),
            'xls.cluster_percent.xls': '{}.eRPKM.xls.cluster_percent.xls'.format(self.prefix),
            'xls.deviation.xls': '{}.eRPKM.xls.deviation.xls'.format(self.prefix),
            'xls.saturation.R': '{}.eRPKM.xls.saturation.R'.format(self.prefix),
            'pdf': '{}.saturation.pdf'.format(self.prefix),
            'r': '{}.saturation.r'.format(self.prefix),
        }

    @toolfuncdeco
    def run(self):
        super(RpkmSaturationTool, self).run()
        self.rpkm_saturation()
        self.rpkm_plot()
        self.set_output()

    @toolfuncdeco
    def rpkm_saturation(self):
        input_file = self.option('bam').path
        output_prefix = self.prefix
        refgene_bed = self.option('bed').path
        percentile_low_bound = self.option('low_bound')
        percentile_up_bound = self.option('up_bound')
        percentile_step = self.option('step')
        rpkm_cutoff = self.option('rpkm_cutof')
        map_qual = self.option('quality')
        cmd = '{} {}'.format(self.program['python'], self.script['RPKM_saturation'])
        cmd += ' -i {}'.format(input_file)
        cmd += ' -o {}'.format(output_prefix)
        cmd += ' -r {}'.format(refgene_bed)
        cmd += ' -l {}'.format(percentile_low_bound)
        cmd += ' -u {}'.format(percentile_up_bound)
        cmd += ' -s {}'.format(percentile_step)
        cmd += ' -c {}'.format(rpkm_cutoff)
        cmd += ' -q {}'.format(map_qual)
        cmd_name = 'run_rpkm_saturation'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        elif command.return_code in [-7]:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    @toolfuncdeco
    def rpkm_plot(self):
        cmd = '{} {}'.format(self.program['perl'], self.script['saturation2plot'])
        cmd += ' -in {}'.format(self.file['xls'])
        cmd += ' -out {}'.format(self.file['xls'])
        cmd_name = 'run_rpkm_plot'
        # runcmd(self, cmd_name, cmd)
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704307")

    @toolfuncdeco
    def set_output(self):
        for source in self.file.values():
            link_name = os.path.join(self.output_dir, os.path.basename(source))
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        else:
            self.end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rpkm_saturation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.mapping.rpkm_saturation',
            'instant': False,
            'options': {
                'bam': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v3/mouse/bam/A2_2.bam',
                'bed': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v3/mouse/GRCm38.96.bed',
                'low_bound': 5,
                'up_bound': 100,
                'step': 5,
                'rpkm_cutof': 0.01,
                'quality': 30
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
