# coding=utf-8
import os
# import unittest
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from decimal import *
from collections import OrderedDict
import copy
__author__ = 'fengyitong'


class ExportCogClassAgent(Agent):
    """
    需输入表达量矩阵及样本分组信息
    """
    def __init__(self, parent):
        super(ExportCogClassAgent, self).__init__(parent)
        options = [
            dict(name="diff_file", type="infile", format="itraq_and_tmt.express_matrix"),
            dict(name="cog_stat", type="infile", format="itraq_and_tmt.common"),
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(ExportCogClassAgent, self).end()


class ExportCogClassTool(Tool):
    def __init__(self, config):
        super(ExportCogClassTool, self).__init__(config)

    def cog_class(self):
        diff_file = self.option('diff_file').prop["path"]
        cmp = os.path.basename(diff_file).split('_diff.xls')[0]
        diff_df = pd.read_csv(diff_file, sep='\t', index_col=0)
        uplist = diff_df[((diff_df['regulate'] == 'up') & (diff_df['significant'] == 'yes'))].index.tolist()
        downlist = diff_df[((diff_df['regulate'] == 'down') & (diff_df['significant'] == 'yes'))].index.tolist()
        cog_df = pd.read_csv(self.option('cog_stat').prop["path"], sep='\t', skiprows=2, usecols=[0, 1, 2, 4], header=None,
                           names=["type", "function_categories", "cog", "cog_list"]).fillna('')
        up_down_cog = pd.DataFrame(
            columns=['Type', 'Functional Categoris', '%s_up_COG' % cmp, '%s_up_COG_list' % cmp,
                     '%s_down_COG' % cmp, '%s_down_COG_list' % cmp])
        de_cog = pd.DataFrame(
            columns=['Type', 'Functional Categoris', '%s_all_COG' % cmp, '%s_all_COG_list' % cmp])
        for index in cog_df.index:
            up_down_dict = dict()
            de_dict = dict()
            cog_accs = cog_df.loc[index, 'cog_list'].split(';')
            up_accs = list(set(uplist) & set(cog_accs))
            down_accs = list(set(downlist) & set(cog_accs))
            if not up_accs and not down_accs:
                continue
            up_down_dict['Type'] = cog_df.loc[index, 'type']
            de_dict['Type'] = cog_df.loc[index, 'type']
            up_down_dict['Functional Categoris'] = cog_df.loc[index, 'function_categories']
            de_dict['Functional Categoris'] = cog_df.loc[index, 'function_categories']
            up_down_dict['%s_up_COG_list' % cmp] = ';'.join(up_accs)
            up_down_dict['%s_up_COG' % cmp] = str(len(up_accs))
            up_down_dict['%s_down_COG_list' % cmp] = ';'.join(down_accs)
            up_down_dict['%s_down_COG' % cmp] = str(len(down_accs))
            de_dict['%s_all_COG_list' % cmp] = ';'.join(up_accs+down_accs)
            de_dict['%s_all_COG' % cmp] = str(len(up_accs+down_accs))
            up_down_cog = up_down_cog.append(up_down_dict,ignore_index=True)
            de_cog = de_cog.append(de_dict, ignore_index=True)
        de_cog_file = os.path.join(self.output_dir, '%s_all_cog_class_table.xls' % cmp)
        up_down_cog_file = os.path.join(self.output_dir, '%s_up_down_cog_class_table.xls' % cmp)
        de_cog.to_csv(de_cog_file, sep='\t', index=False, header=True)
        up_down_cog.to_csv(up_down_cog_file, sep='\t', index=False, header=True)
        return de_cog_file, up_down_cog_file

    def set_output(self):
        # all ready write results to output
        pass

    def run(self):
        super(ExportCogClassTool, self).run()
        self.cog_class()
        self.set_output()
        self.end()



