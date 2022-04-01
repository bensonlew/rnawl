# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SweepGraphAgent(Agent):
    def __init__(self, parent):
        super(SweepGraphAgent, self).__init__(parent)
        options = [
            {"name": "vcftools_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},  # result dir of vcftools_stat
            {'name': 'threshold', 'type': 'float', 'default': 0.05},
        ]
        self.add_option(options)
        self.step.add_steps("sweep_graph")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.sweep_graph.start()
        self.step.update()

    def step_finish(self):
        self.step.sweep_graph.finish()
        self.step.update()

    def check_options(self):
        if not self.option("vcftools_dir").is_set:
            raise OptionError("必须设置输入vcftools结果文件夹")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(SweepGraphAgent, self).end()


class SweepGraphTool(Tool):
    def __init__(self, config):
        super(SweepGraphTool, self).__init__(config)
        self._version = "v1.0"
        self.program = {
            'rscript': 'bioinfo/rna/miniconda2/lib/R/bin/Rscript',
        }
        self.script = {
            'sweep': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/selective_sweep.r'),
        }
        self.file = {
            'outtable': os.path.join(self.output_dir, 'result_table.xls'),
        }

    def run(self):
        super(SweepGraphTool, self).run()
        self.merge_results()
        self.draw_sweep()
        self.set_output()
        self.end()

    def merge_results(self):
        stats_dir = self.option('vcftools_dir').prop['path']
        self.file['fst'] = glob.glob(os.path.join(stats_dir, '*.windowed.weir.fst'))[0]
        groups = os.path.basename(self.file['fst']).split('.')[0]
        group1, group2 = groups.split('_vs_')
        self.file['pi1'] = os.path.join(stats_dir, group1 + '.windowed.pi')
        self.file['pi2'] = os.path.join(stats_dir, group2 + '.windowed.pi')
        self.file['tajima1'] = os.path.join(stats_dir, group1 + '.Tajima.D')
        self.file['tajima2'] = os.path.join(stats_dir, group2 + '.Tajima.D')
        self.file['outgraph'] = os.path.join(self.output_dir, groups)
        # process pi results
        pi1_df = pd.read_table(self.file['pi1'], header=0)
        pi1_df = pi1_df[['CHROM', 'BIN_START', 'PI']]
        pi1_df.rename(columns={'PI': 'pi1'}, inplace=True)
        pi2_df = pd.read_table(self.file['pi2'], header=0)
        pi2_df = pi2_df[['CHROM', 'BIN_START', 'PI']]
        pi2_df.rename(columns={'PI': 'pi2'}, inplace=True)
        pi_df = pi1_df.merge(pi2_df, how='outer', on=['CHROM', 'BIN_START'])
        # process tajima result
        tajima1_df = pd.read_table(self.file['tajima1'], header=0)
        tajima1_df = tajima1_df[['CHROM', 'BIN_START', 'TajimaD']]
        tajima1_df.rename(columns={'TajimaD': 'tajima1'}, inplace=True)
        tajima2_df = pd.read_table(self.file['tajima2'], header=0)
        tajima2_df = tajima2_df[['CHROM', 'BIN_START', 'TajimaD']]
        tajima2_df.rename(columns={'TajimaD': 'tajima2'}, inplace=True)
        tajima_df = tajima1_df.merge(tajima2_df, how='outer', on=['CHROM', 'BIN_START'])
        tajima_df['BIN_START'] = tajima_df.apply(lambda x: x['BIN_START'] + 1, axis=1)
        # process fst result
        fst_df = pd.read_table(self.file['fst'], header=0)
        fst_df = fst_df[['CHROM', 'BIN_START', 'MEAN_FST']]
        fst_df.rename(columns={'MEAN_FST': 'fst'}, inplace=True)
        # merge
        result_table = pi_df.merge(tajima_df, how='outer', on=['CHROM', 'BIN_START'])
        result_table = result_table.merge(fst_df, how='outer', on=['CHROM', 'BIN_START'])
        result_table.rename(columns={'CHROM': 'chr', 'BIN_START': 'pos'},   inplace=True)
        result_table.to_csv(self.file['outtable'], index=False, sep='\t')

    def draw_sweep(self):
        cmd = "{} {} ".format(self.program['rscript'], self.script['sweep'])
        cmd += '--fst {} '.format(self.file['fst'])
        cmd += '--pi1 {} --pi2 {} '.format(self.file['pi1'], self.file['pi2'])
        cmd += '--thre {} --out {} '.format(self.option('threshold'), self.file['outgraph'])
        self.logger.info('start to draw fst+Pi')
        command = self.add_command("sweep_draw", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("绘制selective sweep图完成!")
        else:
            self.set_error("绘制selective sweep图出错！")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        pass