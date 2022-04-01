# -*- coding: utf-8 -*-


from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import json
import pandas as pd
import glob


class TwoSampleDiffModule(Module):
    def __init__(self, work_id):
        super(TwoSampleDiffModule, self).__init__(work_id)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_abun"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_detail", "type": "string"},
            {'name': 'diff_group_name', 'type': 'string', 'default': ''},  # 差异分组
            {'name': 'test_method', 'type': 'string', 'default': 'fisher'},  # 差异检验方法  "chiq", 'fisher'
            {'name': 'side_type', 'type': 'string', 'default': 'two.side'}, # 单尾或双尾检验 two.side,less,greater
            {'name': 'correct', 'type': 'string', 'default':'bonferroni'},
            #创建差异代谢集参数
            {'name': 's2_p_value_fdr', 'type': 'string', 'default': 'P_value'},
            {'name': 's2_p_fdr_condition', 'type': 'string','default': '$gt'},
            {'name': 's2_p_fdr_value', 'type': 'float'},
            {'name': 's2_up_condition', 'type': 'string','default': '$gt'},
            {'name': 's2_up_value', 'type': 'float', 'default': 1.2},
            {'name': 's2_down_condition', 'type': 'string', 'default': '$lt'},
            {'name': 's2_down_value', 'type': 'float', 'default': 0.8},
            {'name': 'mul_metabset', 'type': 'outfile', 'format': 'metabolome.mul_metabset'},
            {'name': 'mix_metabset', "type": "outfile", "format": "metabolome.metabset"}
        ]
        self.add_option(options)
        self.diff_module = self.add_module('metabolome.two_sample_batch')


    def check_options(self):
        if not self.option("metab_table").is_set:
            raise OptionError("请传入丰度矩阵！", code="24700201")
        if not self.option("metab_desc").is_set:
            raise OptionError("请传入metab_desc文件！", code="24700202")
        try:
            self.group_detail = eval(self.option('group_detail'))
        except Exception as e:
            raise OptionError("group_detail  Format Is Wrong")

        self.diff_group = self.option('diff_group_name').split(';')
        if (self.diff_group).sort() != self.group_detail.keys().sort():
            raise OptionError('group detail 和group name 不匹配')

        return True

    def run_two_sample(self):
        self.logger.info("start run diff_sample!")
        options = {
            "metab_table": self.option('metab_table'),
            "group_detail" : self.option("group_detail"),
            "diff_group_name" : self.option('diff_group_name'),
            "test_method" : self.option("test_method"),
            "side_type" : self.option("side_type"),
            "correct" : self.option("correct"),
            'metab_desc': self.option("metab_desc")
        }
        self.diff_module.set_options(options)
        self.diff_module.on('end', self.set_output)
        self.diff_module.run()


    def create_diff_metabset(self,rm_no_name=True):
        if self.option("s2_p_value_fdr") == 'P_value':
            p_col = 'pvalue'
        else:
            p_col = 'corrected_pvalue'
        p_value = self.option('s2_p_fdr_value')

        diff_groups = self.option("diff_group_name").split(';')
        metab_set_map = {}
        for diff in diff_groups:
            files = glob.glob(self.diff_module.output_dir + '/' + diff +'/*_result.xls')
            if files:
                file = files[0]
            else:
                continue

            data = pd.read_table(file, sep='\t',index_col=0)
            if rm_no_name:
                data = data[data['Metabolite']!=data.index]
            if self.option("s2_p_fdr_condition") == '$gt':
                data = data[data[p_col] > p_value]
            elif self.option("s2_p_fdr_condition") == '$gte':
                data = data[data[p_col] >= p_value]

            elif self.option("s2_p_fdr_condition") == '$lte':
                data = data[data[p_col] <= p_value]
            else:
                data = data[data[p_col] < p_value]


            if  self.option("s2_up_condition") == '$gt':
                data1 = data[data['odds_ratio'] > self.option("s2_up_value")]
            elif  self.option("s2_up_condition") == '$gte':
                data1 = data[data['odds_ratio'] >= self.option("s2_up_value")]
            elif self.option("s2_up_condition") == '$lte':
                data1 = data[data['odds_ratio'] <= self.option("s2_up_value")]
            else:
                data1 = data[data['odds_ratio'] < self.option("s2_up_value")]

            if  self.option("s2_down_condition") == '$gt':
                data2 = data[data['odds_ratio'] > self.option("s2_down_value")]
            elif  self.option("s2_down_condition") == '$gte':
                data2 = data[data['odds_ratio'] >= self.option("s2_down_value")]
            elif self.option("s2_down_condition") == '$lte':
                data2 = data[data['odds_ratio'] <= self.option("s2_down_value")]
            else:
                data2 = data[data['odds_ratio'] < self.option("s2_down_value")]

            data_final = pd.concat([data1,data2],axis=0)
            metab_set = list(set(data_final.index.tolist()))
            metab_set_map[diff] = metab_set

        if not os.path.exists(self.work_dir + '/Metabset'):
            os.mkdir(self.work_dir + '/Metabset')

        mix_list = []
        for diff in metab_set_map:
            out_file = self.work_dir + '/Metabset/'+diff+'.metabset.xls'
            with open(out_file, 'w') as fw:
                fw.write('\n'.join(metab_set_map[diff]))
                mix_list.extend(metab_set_map[diff])

        mix_set = set(mix_list)
        if len(mix_set) > 0:
            mul_out = self.work_dir + '/Metabset/mul.metabset.list.xls'
            with open(mul_out, 'w') as fw:
                for diff in metab_set_map:
                    if len(metab_set_map[diff]) > 0:
                        fw.write(diff+'\t'+','.join(metab_set_map[diff])+'\n')

            self.option('mul_metabset', mul_out)

            mix_out = self.work_dir + '/Metabset/Diff_intersection.metabset.xls'
            with open(mix_out,'w') as fw:
                fw.write('metab_id\n')
                fw.write('\n'.join(mix_set)+'\n')
            self.option('mix_metabset', mix_out)


    def set_output(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.diff_module.output_dir, self.output_dir)
        self.create_diff_metabset()
        self.end()

    def run(self):
        super(TwoSampleDiffModule, self).run()
        self.run_two_sample()

    def end(self):
        super(TwoSampleDiffModule, self).end()

