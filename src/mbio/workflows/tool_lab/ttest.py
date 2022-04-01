# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import pandas as pd
from mbio.packages.tool_lab.get_box import get_box_value_row
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
from mbio.packages.tool_lab.common_function import rename_name,rename_name_back
from biocluster.file import download,exists
from biocluster.core.exceptions import OptionError
import HTMLParser
from mbio.packages.tool_lab.common_function import meta_get_file

class TtestWorkflow(Workflow):
    """
    T检验小工具
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TtestWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'tool_lab.table2'},  # 表达矩阵文件
            {"name": "group_file", "type": "infile", "format": "tool_lab.group"},  # 分组文件  meta.otu.group_table
            #{'name': 'group_name', 'type': 'string', 'default': ''},  # 用于检测的组名，eg："A|B"
            {'name': 'group1_name', 'type': 'string', 'default': ''},
            {'name': 'group2_name', 'type': 'string', 'default': ''},
            {'name': 'side_type', 'type': 'string', 'default': 'two-tailed'},
            # 单尾或双尾检验 two-tailed,left-tailed,right-tailed
            {'name': 'mul_test', 'type': 'string', 'default': 'none'},
            # ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "method", "type": "string", "default": "T"}, # 数据是否归一化 T 和 FF
            {"name": "test", "type": "string", "default": "student"}, # 检验方法 student或者welch
            {"name": "ci_method", "type": "string"}, # CI计算方法
            {"name": "ci", "type": "float", "default": 0.95}, # CI值
            {'name': 'project_data', 'type': 'string'},
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")
        if self.option("project_data"):
            self.data_table_raw, self.group_table = meta_get_file(self.option('project_data'), self.work_dir)
        else:
            self.data_table_raw = self.option('table').path
            self.group_table = self.option('group_file').path
        self.file_rename()
        self.run_diff()
        super(TtestWorkflow, self).run()

    def file_rename(self):
        if self.option("method") == "T":
            with open(self.data_table_raw) as f:
                data = f.readlines()
                num = 0
                for i in data:
                    if i.strip():
                        num +=1
                if num == 2:
                    raise OptionError('1行数据无法进行relative数据转换!')
        group_data = pd.read_table(self.group_table, sep='\t')
        group_data.columns = ['#sample', 'group']
        self.group = self.work_dir + '/pick_group.txt'
        groups = [self.option('group1_name'), self.option('group2_name')]
        new_group = group_data[group_data['group'].map(lambda x: x in groups)]
        new_group.to_csv(self.group, sep='\t', header=True, index=False)
        self.data_table = self.work_dir + "/input_table.txt"
        if os.path.exists(self.data_table):
            os.remove(self.data_table)
        self.name_dict = rename_name(self.data_table_raw, self.data_table)

    def run_diff(self):
        #self.this_tool = self.add_tool("tool_lab.diff_test")
        self.this_tool = self.add_tool("tool_lab.metastat2")

        mul_test =  self.option("mul_test")
        if mul_test == 'FDR':
            mul_test = 'fdr'

        side_type_map ={
            "two-tailed" : "two.side",
            "left-tailed" : "less" ,
            "right-tailed" : "greater"
        }
        side_type = side_type_map[self.option("side_type")]

        if self.option("test") == "student":
            options = {
                "student_input": self.data_table,
                "student_group": self.group,
                "student_ci": 0.05,
                "student_correction": mul_test,
                "student_type": side_type,
                "test": self.option("test"),
                "student_gname": 'group',
                "student_coverage": self.option("ci"),
                "student_method": self.option("method"),
            }
        else:
            options = {
                "welch_input": self.data_table,
                "welch_ci": 0.05,
                "welch_group": self.group,
                "welch_correction": mul_test,
                "welch_type": side_type,
                "test": self.option("test"),
                "welch_gname": 'group',
                "welch_coverage": self.option("ci"),
                "welch_method": self.option("method"),
            }
        '''
        options = {
            "student_input": self.option("table").path,
            "student_group": self.group,
            "student_ci": 0.05,
            "student_correction": mul_test,
            "student_type": side_type,
            "test": 'student',
            "student_gname": 'group'
        }
        '''
        self.this_tool.set_options(options)
        self.this_tool.on('end',self.set_db)
        self.this_tool.run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        for file in os.listdir(self.this_tool.output_dir):
            rename_name_back(self.this_tool.output_dir + '/' + file,self.name_dict)
        if self.option("project_data"):
            if os.path.exists(self.this_tool.output_dir + "/input_data"):
                shutil.rmtree(self.this_tool.output_dir + "/input_data")
            os.mkdir(self.this_tool.output_dir + "/input_data")
            os.link(self.data_table_raw, self.this_tool.output_dir + "/input_data/input_table.txt")
            os.link(self.group_table, self.this_tool.output_dir + "/input_data/input_group.txt")
        box_file = self.work_dir + '/box_result.xls'
        get_box_value_row(self.data_table_raw, self.group, box_file,self.option("method"))
        api_name = self.api.api("tool_lab.ttest")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")

        api_name.add_detail_box(main_id,self.this_tool.output_dir,box_file,method=self.option("method"),mul_test=self.option("mul_test"))

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
        super(TtestWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test_sample_diff1',
        'id': 'tsg_3696422',
        'type': 'workflow',
        'options': {
            "table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "group_file": "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/remote_input/group_table/group.txt",
            "group1_name" : "NG_D1",
            "group2_name" : "NG_D3",
            "main_id" : "5e9e6a6017b2bf2049a81be3",

        }
    }

    wsheet = Sheet(data=data)

    wf = TtestWorkflow(wsheet)
    wf.run()