# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class LncrnaLocConsAgent(Agent):
    '''
    last_modify: 2019.04.24
    '''
    def __init__(self, parent):
        super(LncrnaLocConsAgent, self).__init__(parent)
        options = [
            {'name': 'lncrna_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'keep_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'genome', 'type': 'string', 'default': None},
            {'name': 'ways', 'type': 'string', 'default': None},
            {'name': 'thread', 'type': 'int', 'default': 20},
            {'name': 'cutoff', 'type': 'float', 'default': 0.1},
            {'name': 'color', 'type': 'string', 'default': '009933'},
            {'name': 'type_tsv', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'tabular', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('lncrna_loc_cons')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.lncrna_loc_cons.start()
        self.step.update()

    def stepfinish(self):
        self.step.lncrna_loc_cons.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('thread')
        self._memory = '20G'

    def end(self):
        super(LncrnaLocConsAgent, self).end()

class LncrnaLocConsTool(Tool):
    def __init__(self, config):
        super(LncrnaLocConsTool, self).__init__(config)
        self.gtftogenepred = 'bioinfo/align/ucsc_tools/gtfToGenePred'
        self.genepredtobed = 'bioinfo/align/ucsc_tools/genePredToBed'
        self.python = 'program/Python/bin/python'
        self.bwtool = 'bioinfo/lnc_rna/bwtool/bwtool-master/build/bin/bwtool'
        self.rscript = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/miniconda2/bin/Rscript')
        self.filter_bed_by_name_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_bed_by_name.py')
        self.change_bed_chrom_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/change_bed_chrom.py')
        self.plot_loc_cons_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/plot_loc_cons.py')
        self.tsv2line_r = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/tsv2line.r')
        self.loc_cons_statistics_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/loc_cons_statistics.py')
        self.lncrna_genepred = os.path.join(self.work_dir, 'lncrna.genepred')
        self.lncrna_bed = os.path.join(self.work_dir, 'lncrna.bed')
        self.map_txt = os.path.join(
            self.config.SOFTWARE_DIR, 'database/lnc_rna/phastCons/map/{}.txt'.format(self.option('genome'))
        )
        self.filter_bed = os.path.join(self.work_dir, 'filter.bed')
        self.regions_bed = os.path.join(self.work_dir, '{}.bed'.format(self.option('genome')))
        self.phast_cons_bw = os.path.join(self.config.SOFTWARE_DIR,
            'database/lnc_rna/phastCons/{}.phastCons{}way.bw'.format(self.option('genome'), self.option('ways'))
        )
        self.extract_txt = os.path.join(self.work_dir, 'extract.txt')
        self.retain_list = os.path.join(self.work_dir, 'retain.list')
        self.summary_txt = os.path.join(self.work_dir, 'summary.txt')
        self.tabular = os.path.join(self.output_dir, 'lncRNA_location_conservation.tabular')

    def run(self):
        super(LncrnaLocConsTool, self).run()
        self.run_gtftogenepred()
        self.run_genepredtobed()
        self.run_filter_bed_by_name()
        self.run_change_bed_chrom()
        self.run_bwtool_extract()
        self.run_plot_loc_cons()
        self.run_bwtool_summary()
        self.run_loc_cons_statistics()
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
        cmd += ' -o {}'.format(self.regions_bed)
        cmd_name = 'run_change_bed_chrom'
        self.run_code(cmd_name, cmd)

    def run_bwtool_extract(self):
        cmd = '{} extract bed {} {} {}'.format(self.bwtool, self.regions_bed, self.phast_cons_bw, self.extract_txt)
        cmd_name = 'run_bwtool_extract'
        self.run_code(cmd_name, cmd)

    def run_plot_loc_cons(self):
        cmd = '{} {}'.format(self.python, self.plot_loc_cons_py)
        cmd += ' -i {}'.format(self.extract_txt)
        cmd += ' -t {}'.format(self.option('thread'))
        cmd += ' -n {}'.format(self.option('cutoff'))
        cmd += ' -r {}'.format(self.rscript)
        cmd += ' -s {}'.format(self.tsv2line_r)
        cmd += ' -c {}'.format(self.option('color'))
        cmd += ' -k {}'.format(self.retain_list)
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_plot_loc_cons'
        self.run_code(cmd_name, cmd)

    def run_bwtool_summary(self):
        cmd = '{} summary -with-sum -keep-bed -header {} {} {}'.format(
            self.bwtool, self.regions_bed, self.phast_cons_bw, self.summary_txt
        )
        cmd_name = 'run_bwtool_summary'
        self.run_code(cmd_name, cmd)

    def run_loc_cons_statistics(self):
        cmd = '{} {}'.format(self.python, self.loc_cons_statistics_py)
        cmd += ' -i {}'.format(self.summary_txt)
        cmd += ' -r {}'.format(self.retain_list)
        cmd += ' -t {}'.format(self.option('type_tsv').path)
        cmd += ' -o {}'.format(self.tabular)
        cmd_name = 'run_loc_cons_statistics'
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
            'id': 'lncrna_loc_cons_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.lncrna_loc_cons',
            'instant': False,
            'options': {
                'lncrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33915/upload_results/Filtered_result/filtered_file/all_lncrna.gtf',
                'keep_list': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33915/upload_results/Filtered_result/filtered_file/known_lncrna_ids.list',
                'genome': 'hg38',
                'ways': '4',
                'type_tsv': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33915/upload_results/Filtered_result/filtered_file/trans_type.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
