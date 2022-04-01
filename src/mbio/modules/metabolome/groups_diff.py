# -*- coding: utf-8 -*-


from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import json
import pandas as pd
import glob


class GroupsDiffModule(Module):
    def __init__(self, work_id):
        super(GroupsDiffModule, self).__init__(work_id)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_abun"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_name", "type": "string",'default':"group"},
            {"name": "test_method", "type" : "string", "default": "ow"},
            {'name': 'post_hoc', 'type': 'string', 'default': 'scheffe'},
            {'name': 'coverage', 'type': 'float', 'default':0.95}
        ]
        self.add_option(options)
        self.diff_tool = self.add_tool("metabolome.diff.metastat")

    def check_options(self):
        if not self.option("metab_table").is_set:
            raise OptionError("请传入丰度矩阵！", code="24700201")
        if not self.option("group").is_set:
            raise OptionError("请传入group！", code="24700201")
        if not self.option("metab_desc").is_set:
            raise OptionError("请传入metab_desc！", code="24700201")

        return True

    def run_groups_fun(self):
        exp_file = self.option('metab_table').path
        if self.option("test_method") == "ow":
            options = {
                "anova_input": exp_file,
                "anova_group": self.group,
                "anova_correction": 'fdr',
                "test": 'anova',
                "anova_gname": self.option("group_name"),
                "anova_methor": self.option("post_hoc"),
                "anova_coverage": self.option("coverage")
            }
        elif self.option("test_method") == "kw":
            options = {
                "kru_H_input": exp_file,
                "kru_H_group": self.group,
                "kru_H_correction": "fdr",
                "test": 'kru_H',
                "kru_H_gname": self.option("group_name"),
                "kru_H_methor": self.option("post_hoc"),
                "kru_H_coverage": self.option("coverage")
            }
        self.diff_tool.set_options(options)
        self.diff_tool.on('end',self.set_output)
        self.diff_tool.run()


    def run(self):
        super(GroupsDiffModule, self).run()
        ##修改group file的表头
        group_data = pd.read_table(self.option("group").path,sep='\t')
        group_data.columns = ['#sample',self.option('group_name')]
        self.group = self.work_dir+'/new_header_group.txt'
        group_data.to_csv(self.group,sep='\t',header=True,index=False)
        self.run_groups_fun()



    def set_output(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.diff_tool.output_dir , self.output_dir)
        target_dir = self.output_dir
        metab_desc = self.option("metab_desc").path
        self.add_metab_name_to_result_file(target_dir, metab_desc)
        self.end()

    def add_metab_name_to_result_file(self,tar_dir,desc_file):
        desc = pd.read_table(desc_file,sep='\t',index_col=0)
        metab_name = desc['Metabolite']
        results = glob.glob(tar_dir+'/*_result.xls')
        for r_file in results:
            data = pd.read_table(r_file,sep='\t',index_col=0)
            data = data.join(metab_name)
            data.to_csv(r_file,sep='\t')

    def end(self):
        super(GroupsDiffModule, self).end()

