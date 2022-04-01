# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.tool_lab.get_box import get_box_value_row
import os
from mbio.packages.tool_lab.common_function import rename_name,rename_name_back
from biocluster.file import download,exists
import HTMLParser
from mbio.packages.tool_lab.common_function import meta_get_file
import shutil


class AnovaWorkflow(Workflow):
    """
    单因素方差小工具
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnovaWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'tool_lab.table2'},  # 表达矩阵文件
            {"name": "group_file", "type": "infile", "format": "tool_lab.simple"},  # 分组文件
            {'name': 'group_name', 'type': 'string', 'default': ''},  # 用于检测的组名，eg："A|B|C"  或 A,B,C
            {'name': 'mul_test', 'type': 'string', 'default': 'none'},
            # ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "method", "type": "string", "default": "none"},  # 数据转换
            {"name": "methor", "type": "string", "default": "tukeykramer"},  # Post-hoc检验方法
            {"name": "coverage", "type": "float", "default": 0.95},  # Post-hoc检验值
            {'name': 'project_data', 'type': 'string'},
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_file(self):
        if self.option("method") == "T":
            with open(self.data_table_raw) as f:
                data = f.readlines()
                num = 0
                for i in data:
                    if i.strip():
                        num +=1
                if num == 2:
                    raise OptionError('1行数据无法进行relative数据转换!')
        if '|' in self.option('group_name'):
            self.groups = self.option('group_name').split('|')
        elif ',' in self.option('group_name'):
            self.groups = self.option('group_name').split(',')
        elif '，' in self.option('group_name'):
            self.groups = self.option('group_name').split('，')
        else:
            self.set_error("组名之间要用，或|分割")
        if len(self.groups) < 3:
            self.set_error("分组数小于3组")

        self.groups = [g.strip() for g in self.groups]


    def group_check(self):
        if os.path.getsize(self.group_table) == 0:
            self.set_error('分组文件为空')

        if os.path.getsize(self.data_table_raw) == 0:
            self.set_error('数据表文件为空')

        group_data = pd.read_table(self.group_table, sep='\t',header=0)
        table =  pd.read_table(self.data_table_raw, sep='\t',header=0)
        if len(group_data)==0 :
            self.set_error("分组文件为空")
        if len(table)==0 :
            self.set_error("table文件为空")

        for sample in group_data[group_data.columns[0]]:
            if sample not in table.columns:
                self.set_error('分组文件的样本%s 不在另一文件中'%sample)
        for g in self.groups:
            if g not in group_data[group_data.columns[1]].tolist():
                self.set_error('输入的组名没有在分组文件中')

        self.data_table = self.work_dir + "/input_table.txt"
        if os.path.exists(self.data_table):
            os.remove(self.data_table)
        self.name_dict = rename_name(self.data_table_raw, self.data_table)


    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")
        if self.option("project_data"):
            self.data_table_raw, self.group_table = meta_get_file(self.option('project_data'), self.work_dir)
        else:
            self.data_table_raw = self.option('table').path
            self.group_table = self.option('group_file').path
        self.check_file()
        self.group_check()
        self.run_diff()

        super(AnovaWorkflow, self).run()

    def run_diff(self):
        group_data = pd.read_table(self.group_table, sep='\t')
        group_data.columns = ['#sample', 'group']
        self.group = self.work_dir + '/pick_group.txt'
        new_group = group_data[group_data['group'].map(lambda x: x in self.groups)]
        new_group.to_csv(self.group, sep='\t', header=True, index=False)

        self.this_tool = self.add_tool("tool_lab.metastat2")

        mul_test = self.option("mul_test")
        if mul_test == 'FDR':
            mul_test = 'fdr'

        options = {
                "anova_input": self.data_table,
                "anova_group": self.group,
                "anova_correction": mul_test,
                "test": 'anova',
                "anova_gname": 'group',
                "anova_methor": self.option("methor"),
                "anova_coverage": self.option("coverage"),
                "anova_method": self.option("method"),
            }
        self.this_tool.set_options(options)
        self.this_tool.on('end',self.set_db)
        self.this_tool.run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        if os.path.exists(self.this_tool.output_dir + '/anova_result.xls'):
            os.rename(self.this_tool.output_dir + '/anova_result.xls',self.this_tool.output_dir + '/anova_result_bak1.xls')
            with open(self.this_tool.output_dir + '/anova_result_bak1.xls') as f,open(self.this_tool.output_dir + '/anova_result.xls',"w") as t:
                data = f.readlines()
                if data[0].strip().split("\t")[0] != "name":
                    t.write("name"+"\t"+data[0].strip()+"\n")
                else:
                    t.write(data[0])
                for i in data[1:]:
                    t.write(i)
                os.remove(self.this_tool.output_dir + '/anova_result_bak1.xls')
        if os.path.exists(self.this_tool.output_dir + '/anova_boxfile.xls'):
            os.rename(self.this_tool.output_dir + '/anova_boxfile.xls',self.this_tool.output_dir + '/anova_boxfile_bak1.xls')
            with open(self.this_tool.output_dir + '/anova_boxfile_bak1.xls') as f,open(self.this_tool.output_dir + '/anova_boxfile.xls',"w") as t:
                data = f.readlines()
                if data[0].strip().split("\t")[0] != "name":
                    t.write("name"+"\t"+data[0].strip()+"\n")
                else:
                    t.write(data[0])
                for i in data[1:]:
                    t.write(i)
                os.remove(self.this_tool.output_dir + '/anova_boxfile_bak1.xls')
        for file in os.listdir(self.this_tool.output_dir):
            rename_name_back(self.this_tool.output_dir + '/' + file,self.name_dict)
        box_file = self.work_dir + '/box_result.xls'
        get_box_value_row(self.data_table_raw, self.group, box_file,self.option('method'))
        api_name = self.api.api("tool_lab.anova")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        if self.option("group_name").find(";") >= 1:
            group_names = self.option("group_name").split(';')
        elif self.option("group_name").find("|") >= 1:
            group_names = self.option("group_name").split('|')
        else:
            group_names = self.option("group_name").split(',')
        api_name.add_detail_box(main_id,self.this_tool.output_dir,box_file,self.option('methor'),self.option('method'),self.option('mul_test'),group_names=group_names)
        if self.option("project_data"):
            if os.path.exists(self.this_tool.output_dir + "/input_data"):
                shutil.rmtree(self.this_tool.output_dir + "/input_data")
            os.mkdir(self.this_tool.output_dir + "/input_data")
            os.link(self.data_table_raw, self.this_tool.output_dir + "/input_data/input_table.txt")
            os.link(self.group_table, self.this_tool.output_dir + "/input_data/input_group.txt")
        self.end()


    def end(self):
        result_dir = self.add_upload_dir(self.this_tool.output_dir)
        relpath_rules = [
            [".", "", "比较结果文件夹", 0, ],
        ]
        regexps = [
            [r".*/.*_result.xls", "xls", "显著性比较结果表，包括均值，标准差，p值", 0],
            [r".*/.*_CI.xls", "xls", "两样本比较的置信区间值", 0]
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(AnovaWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test_sample_diff',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "group_file": "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/remote_input/group_table/group.txt",
            "group_name" : "NG_D1|NG_D3|SMG_D1",
            "main_id" : "5e9e6a6017b2bf2049a81be3",

        }
    }

    wsheet = Sheet(data=data)
    wf = AnovaWorkflow(wsheet)
    wf.run()