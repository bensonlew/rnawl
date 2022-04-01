# -*- coding: utf-8 -*-


from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import json
import pandas as pd


class TwoSampleBatchModule(Module):
    def __init__(self, work_id):
        super(TwoSampleBatchModule, self).__init__(work_id)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_abun"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_detail", "type": "string"},
            {'name': 'diff_group_name', 'type': 'string', 'default': ''},  # 差异分组
            {'name': 'test_method', 'type': 'string', 'default': 'fisher'},  # 差异检验方法  "chiq", 'fisher'
            {'name': 'side_type', 'type': 'string', 'default': 'two.side'}, # 单尾或双尾检验 two.side,less,greater
            {'name': 'correct', 'type': 'string', 'default':'bonferroni'}
        ]
        self.add_option(options)

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

    def run_two_sample(self,sample1, sample2):
        two_sample_tool = self.add_tool("metabolome.diff.metastat")
        if self.option("test_method") == "chiq":
            options = {
                "chi_input": self.option("metab_table").path,
                "chi_sample1": sample1,
                "chi_sample2": sample2,
                "chi_correction": self.option("correct"),
                "test": 'chi',
                "chi_methor": 'DiffBetweenPropAsymptotic',
                "chi_coverage": 0.95,
            }
        else:
            options = {
                "fisher_input": self.option("metab_table").path,
                "fisher_ci": 0.05,
                "fisher_sample1": sample1,
                "fisher_sample2": sample2,
                "fisher_correction": self.option("correct"),
                "test": 'fisher',
                "fisher_type": self.option("side_type"),
                "fisher_methor": "DiffBetweenPropAsymptotic",
                "fisher_coverage": 0.95
            }
        two_sample_tool.set_options(options)
        return two_sample_tool

    def run_two_sample_batch(self):
        for diff_group in self.diff_group:
            groups = diff_group.split('_vs_')
            g1 = groups[0]
            g2 = groups[1]
            s1 = self.group_detail[g1][0]
            s2 = self.group_detail[g2][0]
            two_sample_tool = self.run_two_sample(s1,s2)
            self.tool_list.append(two_sample_tool)
            self.diff_group_map_tool[diff_group] = two_sample_tool

        self.on_rely(self.tool_list, self.set_output)

        with open(self.work_dir+'/tool_map.txt','w') as fw:
            for diff_group in self.diff_group_map_tool:
                fw.write(diff_group+'\t'+self.diff_group_map_tool[diff_group]._id+'\n')

        for tool in self.tool_list:
            tool.run()

    def set_output(self):

        desc = pd.read_table(self.option("metab_desc").path, sep='\t',index_col=0)
        metab_names = pd.DataFrame(desc['Metabolite'])
        for diff_group in self.diff_group_map_tool:
            target_dir =  self.output_dir + '/' +  diff_group
            current_tool = self.diff_group_map_tool[diff_group]
            if os.path.exists(target_dir):
                shutil.rmtree(target_dir)
            shutil.copytree(current_tool.output_dir, target_dir)
            ## 计算比值
            if  self.option("test_method") == "chiq":
                result = target_dir+'/chi_result.xls'
            else:
                result = target_dir+'/fisher_result.xls'
            data = pd.read_table(result,sep='\t',index_col=0)
            col1 = data.columns[0]
            col2 = data.columns[1]

            data['odds_ratio'] = data[col1]/data[col2]
            mv = data['odds_ratio'][data['odds_ratio']!=float('inf')].max()
            data['odds_ratio'] = data['odds_ratio'].replace(float('inf'), mv*2)

            data.to_csv(result,sep='\t')
            ## add metab name
            for file_name in os.listdir(target_dir):
                file = target_dir + '/' + file_name
                data = pd.read_table(file, sep='\t',index_col=0)
                #out = pd.concat([data, metab_names],axis=0)
                #or
                #pd.merge(data.reset_index(),metab_names,left_on='index',right_on='metab_id')
                out = data.join(metab_names)
                out.to_csv(file, sep='\t')

        self.end()

    def run(self):
        super(TwoSampleBatchModule, self).run()
        self.tool_list = []
        self.diff_group_map_tool = {}
        self.run_two_sample_batch()

    def end(self):
        super(TwoSampleBatchModule, self).end()

