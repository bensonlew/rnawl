# coding=utf-8
import os
# import unittest
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from decimal import *
from collections import OrderedDict
__author__ = 'fengyitong'


class ExportGoClassAgent(Agent):
    """
    需输入表达量矩阵及样本分组信息
    """
    def __init__(self, parent):
        super(ExportGoClassAgent, self).__init__(parent)
        options = [
            dict(name="diff_file", type="infile", format="itraq_and_tmt.express_matrix"),
            dict(name="go_stat", type="infile", format="itraq_and_tmt.common"),
            dict(name="go_version", type="string", default="2019"),

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(ExportGoClassAgent, self).end()


class ExportGoClassTool(Tool):
    def __init__(self, config):
        super(ExportGoClassTool, self).__init__(config)

    def go_class(self):
        diff_file = self.option('diff_file').prop["path"]
        cmp = os.path.basename(diff_file).split('_diff.xls')[0]
        diff_df = pd.read_csv(diff_file, sep='\t', index_col=0)
        uplist = diff_df[((diff_df['regulate'] == 'up') & (diff_df['significant'] == 'yes'))].index.tolist()
        downlist = diff_df[((diff_df['regulate'] == 'down') & (diff_df['significant'] == 'yes'))].index.tolist()
        go_df = pd.read_csv(self.option('go_stat').prop["path"], sep='\t').fillna('')
        de_go = pd.DataFrame(columns=['Term type','Term','GO','%s_all num'%cmp,'%s_all percent'%cmp,'%s_all list' % cmp])
        up_down_go = pd.DataFrame(columns=['Term type','Term','GO','%s_up num' % cmp,'%s_up percent' % cmp,'%s_up list' % cmp,'%s_down num' % cmp,
                                           '%s_down percent' % cmp,'%s_down list' % cmp])
        for term in ["biological_process", "cellular_component", "molecular_function"]:
            tmp_go = go_df[go_df['GO (Lev1)'] == term]
            # for n in range(tmp_go.shape[0]):
            for index in tmp_go.index:
                de_dict = OrderedDict()
                up_down_dict = OrderedDict()
                go_accs = tmp_go.loc[index, 'Seq List'].split(';')
                if not go_accs:
                    continue
                up_accs = list(set(uplist) & set(go_accs))
                down_accs = list(set(downlist) & set(go_accs))
                if not up_accs and not down_accs:
                    continue
                #避免分母为0,by xuxi 20211014
                if len(uplist + downlist) == 0:
                    len_uplist_downlist = 9999
                else:
                    len_uplist_downlist = len(uplist + downlist)
                if len(uplist) == 0:
                    len_uplist = 9999
                else:
                    len_uplist = len(uplist)
                if len(downlist) == 0:
                    len_downlist = 9999
                else:
                    len_downlist = len(downlist)
                #
                de_dict['Term type'] = term
                up_down_dict['Term type'] = term
                de_dict['Term'] = tmp_go.loc[index, 'GO Term (Lev2)']
                up_down_dict['Term'] = tmp_go.loc[index, 'GO Term (Lev2)']
                de_dict['GO'] = tmp_go.loc[index, 'GO ID (Lev2)']
                up_down_dict['GO'] = tmp_go.loc[index, 'GO ID (Lev2)']
                de_dict['%s_all num'%cmp] = len(up_accs + down_accs)
                de_dict['%s_all percent'%cmp] = str(float(len(up_accs + down_accs)) / len_uplist_downlist) + '(' + str(len(up_accs + down_accs)) + '/' + str(len(uplist + downlist)) + ')'
                de_dict['%s_all list' % cmp] = ';'.join(up_accs + down_accs)
                up_down_dict['%s_up num' % cmp] = len(up_accs)
                # up_down_dict['%s_up percent' % cmp] = float(len(up_accs)) / len(go_accs)
                up_down_dict['%s_up percent' % cmp] = str(float(len(up_accs)) / len_uplist) + '(' + str(len(up_accs)) + '/' + str(len(uplist)) + ')'
                up_down_dict['%s_up list' % cmp] = ';'.join(up_accs)
                up_down_dict['%s_down num' % cmp] = len(down_accs)
                # up_down_dict['%s_down percent' % cmp] = float(len(down_accs)) / len(go_accs)
                up_down_dict['%s_down percent' % cmp] = str(float(len(down_accs)) / len_downlist) + '(' + str(len(down_accs)) + '/' + str(len(downlist)) + ')'
                up_down_dict['%s_down list' % cmp] = ';'.join(down_accs)
                de_go = de_go.append(de_dict,ignore_index=True)
                up_down_go = up_down_go.append(up_down_dict,ignore_index=True)
        de_go_file = os.path.join(self.output_dir,'%s_all_go_class_table.xls'%cmp)
        up_down_go_file = os.path.join(self.output_dir,'%s_up_down_go_class_table.xls'%cmp)
        de_go.to_csv(de_go_file,sep='\t', index=False, header=True)
        up_down_go.to_csv(up_down_go_file,sep='\t', index=False, header=True)
        return de_go_file, up_down_go_file

    def set_output(self):
        # all ready write results to output
        pass

    def run(self):
        super(ExportGoClassTool, self).run()
        self.go_class()
        self.set_output()
        self.end()



