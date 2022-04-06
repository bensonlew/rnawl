# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class LncrnaSeqConsAgent(Agent):
    '''
    last_modify: 2019.04.02
    '''
    def __init__(self, parent):
        super(LncrnaSeqConsAgent, self).__init__(parent)
        options = [
            {'name': 'lncrna_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'keep_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'old_genome', 'type': 'string', 'default': None},
            {'name': 'new_genome', 'type': 'string', 'default': None},
            {'name': 'type_tsv', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'tabular', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('lncrna_seq_cons')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.lncrna_seq_cons.start()
        self.step.update()

    def stepfinish(self):
        self.step.lncrna_seq_cons.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '8G'

    def end(self):
        super(LncrnaSeqConsAgent, self).end()

class LncrnaSeqConsTool(Tool):
    def __init__(self, config):
        super(LncrnaSeqConsTool, self).__init__(config)
        self.gtftogenepred = 'bioinfo/align/ucsc_tools/gtfToGenePred'
        self.genepredtobed = 'bioinfo/align/ucsc_tools/genePredToBed'
        self.python = 'miniconda2/bin/python'
        self.liftover = 'bioinfo/align/ucsc_tools/liftOver'
        self.bedtools = 'bioinfo/rna/bedtools2-master/bin/bedtools'
        self.filter_bed_by_name_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_bed_by_name.py')
        self.change_bed_chrom_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/change_bed_chrom.py')
        self.seq_cons_statistics_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/seq_cons_statistics.py')
        self.lncrna_genepred = os.path.join(self.work_dir, 'lncrna.genepred')
        self.lncrna_bed = os.path.join(self.work_dir, 'lncrna.bed')
        self.filter_bed = os.path.join(self.work_dir, 'filter.bed')
        self.map_txt = os.path.join(
            self.config.SOFTWARE_DIR, 'database/lnc_rna/liftOver/map/{}.txt'.format(self.option('old_genome'))
        )
        self.old_bed = os.path.join(self.work_dir, '{}.bed'.format(self.option('old_genome')))
        self.map_chain = os.path.join(
            self.config.SOFTWARE_DIR, 'database/lnc_rna/liftOver/{}To{}.over.chain.gz'.format(
                self.option('old_genome'), self.option('new_genome').capitalize()
            )
        )
        self.new_bed = os.path.join(self.work_dir, '{}.bed'.format(self.option('new_genome')))
        self.unmapped = os.path.join(self.work_dir, 'unMapped')
        self.ref_bed = os.path.join(
            self.config.SOFTWARE_DIR, 'database/lnc_rna/liftOver/bed/{}.bed'.format(self.option('new_genome'))
        )
        self.intersect_tsv = os.path.join(self.work_dir, 'intersect.tsv')
        self.tabular = os.path.join(self.output_dir, 'lncRNA_sequence_conservation.tabular')

    def run(self):
        super(LncrnaSeqConsTool, self).run()
        self.run_gtftogenepred()
        self.run_genepredtobed()
        self.run_filter_bed_by_name()
        self.run_change_bed_chrom()
        self.run_liftover()
        self.run_bedtools()
        self.run_seq_cons_statistics()
        self.set_output()
        self.end()

    def run_gtftogenepred(self):
        cmd = '{} {} {}'.format(self.gtftogenepred, self.option('lncrna_gtf').path, self.lncrna_genepred)
        cmd_name = 'run_gtftogenepred'
        self.run_code(cmd_name, cmd)

    def run_genepredtobed(self):
        cmd = '{} {} {}'.format(self.genepredtobed, self.lncrna_genepred, self.lncrna_bed)
        cmd_name = 'run_genepredtobed'
        self.run_code(cmd_name, cmd)

    def run_filter_bed_by_name(self):
        cmd = '{} {}'.format(self.python, self.filter_bed_by_name_py)
        cmd += ' -i {}'.format(self.lncrna_bed)
        cmd += ' -n {}'.format(self.option('keep_list').path)
        cmd += ' -o {}'.format(self.filter_bed)
        cmd_name = 'run_filter_bed_by_name'
        self.run_code(cmd_name, cmd)

    def run_change_bed_chrom(self):
        cmd = '{} {}'.format(self.python, self.change_bed_chrom_py)
        cmd += ' -i {}'.format(self.filter_bed)
        cmd += ' -m {}'.format(self.map_txt)
        cmd += ' -o {}'.format(self.old_bed)
        cmd_name = 'run_change_bed_chrom'
        self.run_code(cmd_name, cmd)

    def run_liftover(self):
        cmd = '{} {} {} {} {}'.format(self.liftover, self.old_bed, self.map_chain, self.new_bed, self.unmapped)
        cmd_name = 'run_liftover'
        self.run_code(cmd_name, cmd)

    def run_bedtools(self):
        cmd = '{} intersect -s -wo'.format(self.bedtools)
        cmd += ' -a {}'.format(self.new_bed)
        cmd += ' -b {}'.format(self.ref_bed)
        cmd += ' > {}'.format(self.intersect_tsv)
        cmd_name = 'run_bedtools'
        self.run_code(cmd_name, cmd, shell=True)

    def run_seq_cons_statistics(self):
        cmd = '{} {}'.format(self.python, self.seq_cons_statistics_py)
        cmd += ' -l {}'.format(self.new_bed)
        cmd += ' -b {}'.format(self.intersect_tsv)
        cmd += ' -t {}'.format(self.option('type_tsv').path)
        cmd += ' -o {}'.format(self.tabular)
        cmd_name = 'run_seq_cons_statistics'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run {}, abord'.format(n))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('tabular').set_path(self.tabular)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'lncrna_seq_cons_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.lncrna_seq_cons',
            'instant': False,
            'options': {
                'lncrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33915/upload_results/Filtered_result/filtered_file/all_lncrna.gtf',
                'keep_list': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33915/upload_results/Filtered_result/filtered_file/known_lncrna_ids.list',
                'old_genome': 'hg38',
                'new_genome': 'mm10',
                'type_tsv': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33915/upload_results/Filtered_result/filtered_file/trans_type.xls'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
