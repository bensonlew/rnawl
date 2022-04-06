# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
from biocluster.tool import Tool
import shutil
import unittest

class GffcompareAgent(Agent):
    '''
    last_modify: 2019.01.24
    '''
    def __init__(self, parent):
        super(GffcompareAgent, self).__init__(parent)
        options = [
            {'name': 'merged_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'gffcmp_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'gffcmp_tmap', 'type': 'outfile', 'format': 'lnc_rna.tmap'},
        ]
        self.add_option(options)
        self.step.add_steps('gffcompare')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.gffcompare.start()
        self.step.update()

    def stepfinish(self):
        self.step.gffcompare.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('merged_gtf', self.option('merged_gtf').prop['path']))
        if not self.option('merged_gtf').is_set:
            raise OptionError('merged GTF file must be provided')
        self.logger.debug('{} - {}'.format('ref_gtf', self.option('ref_gtf').prop['path']))
        if not self.option('ref_gtf').is_set:
            raise OptionError('reference annotation GTF must be provided')
        self.logger.debug('{} - {}'.format('gffcmp_gtf', self.option('gffcmp_gtf')))
        self.logger.debug('{} - {}'.format('gffcmp_tmap', self.option('gffcmp_tmap')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('merged_gtf').prop['path']) / 1024 ** 3 + 4)

    def end(self):
        super(GffcompareAgent, self).end()

class GffcompareTool(Tool):
    def __init__(self, config):
        super(GffcompareTool, self).__init__(config)
        self.gffcompare = 'bioinfo/rna/gffcompare-0.9.8.Linux_x86_64/gffcompare'
        self.python = 'miniconda2/bin/python'
        self.filter_gtf_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_gtf.py')
        self.merged_gtf = os.path.join(self.work_dir, 'merged.gtf')
        self.outprefix = os.path.join(self.work_dir, 'gffcmp')
        self.gffcmp_annotated_gtf = os.path.join(self.work_dir, 'gffcmp.annotated.gtf')
        self.gffcmp_filter_gtf = os.path.join(self.work_dir, 'gffcmp.filter.gtf')
        self.gffcmp_merged_gtf_tmap = os.path.join(self.work_dir, 'gffcmp.merged.gtf.tmap')

    def run(self):
        super(GffcompareTool, self).run()
        self.run_gffcompare()
        self.run_filter_gtf()
        self.set_output()
        self.end()

    def run_gffcompare(self):
        shutil.copy(self.option('merged_gtf').prop['path'], self.merged_gtf)
        cmd = '{} {} -V '.format(self.gffcompare, self.merged_gtf)
        cmd += '-r {} '.format(self.option('ref_gtf').prop['path'])
        cmd += '-o {}'.format(self.outprefix)
        cmd_name = 'run_gffcompare'
        self.run_code(cmd_name, cmd)

    def run_filter_gtf(self):
        cmd = '{} {} '.format(self.python, self.filter_gtf_py)
        cmd += '--ref {} '.format(self.option('ref_gtf').prop['path'])
        cmd += '--raw {} '.format(self.gffcmp_annotated_gtf)
        cmd += '--filter {}'.format(self.gffcmp_filter_gtf)
        cmd_name = 'run_filter_gtf'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
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

        gffcmp_gtf = os.path.join(self.output_dir, 'gffcmp.gtf')
        if os.path.exists(gffcmp_gtf):
            os.remove(gffcmp_gtf)
        os.link(self.gffcmp_filter_gtf, gffcmp_gtf)
        self.logger.info('succeed in linking {} to {}'.format(self.gffcmp_filter_gtf, gffcmp_gtf))
        self.option('gffcmp_gtf', gffcmp_gtf)

        gffcmp_tmap = os.path.join(self.output_dir, 'gffcmp.tmap')
        if os.path.exists(gffcmp_tmap):
            os.remove(gffcmp_tmap)
        os.link(self.gffcmp_merged_gtf_tmap, gffcmp_tmap)
        self.logger.info('succeed in linking {} to {}'.format(self.gffcmp_merged_gtf_tmap, gffcmp_tmap))
        self.option('gffcmp_tmap', gffcmp_tmap)

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
            'id': 'gffcompare_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.assemble.gffcompare',
            'instant': False,
            'options': {
                'merged_gtf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/assemble/merged.gtf',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()