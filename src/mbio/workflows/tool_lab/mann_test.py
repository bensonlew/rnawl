# !usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import pandas as pd
import os
import subprocess
import re
from mbio.packages.tool_lab.common_function import rename_name,rename_name_back
from biocluster.file import download,exists
import HTMLParser
from mbio.packages.tool_lab.common_function import meta_get_file
import shutil

class MannTestWorkflow(Workflow):
    """
    秩和检验小工具
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MannTestWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "mann_input", "type": "infile", "format": "tool_lab.table2"},  # 秩和检验的输入文件
            {"name": "mann_group", "type": "infile", "format": "tool_lab.group_table"},  # 秩和检验的输入分组文件
            {"name": "sample", "type": "string", "default": "column"},  # 样本名为列标签或行标签
            {"name": "compare_type", "type": "string", "default": "multi"},  # 比较策略
            {"name": "group_name1", "type": "string"},  # 秩和检验样本组1名称
            {"name": "group_name2", "type": "string"},  # 秩和检验样本组2名称
            {"name": "group_name", "type": "string"},  # 秩和检验样本组名称
            {"name": "mann_ci", "type": "float", "default": 0.05},  # 秩和检验的显著性水平
            {"name": "correction", "type": "string", "default": "none"},  # 秩和检验的多重检验校正
            {"name": "mann_type", "type": "string", "default": "double"},  # 秩和检验的选择单尾或双尾检验
            {"name": "kru_methor", "type": "string", "default": "tukeykramer"},  # 多组检验方法
            {"name": "coverage", "type": "float", "default": 0.95},  # 秩和检验的置信度
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "method", "type": "string", "default": "T"},  # 数据是否归一化 无：FF ，relative：T
            {'name': 'project_data', 'type': 'string'},
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.mann = self.add_tool("tool_lab.mann_test")

    def check_options(self):
        """
        检查参数
        """
        #if not self.option("mann_input").is_set:
        #    raise OptionError('必须设置秩和检验输入的文件')
        #if not self.option("mann_group").is_set:
        #    raise OptionError('必须设置秩和检验输入的分组文件')
        if self.option("correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
            raise OptionError('该多重检验校正的方法不被支持')
        if self.option("kru_methor") not in ["tukeykramer", "gameshowell", "welchuncorrected", "scheffe","dunn", "nemenyi", "conover-iman", "steel_dwass"]:
            raise OptionError('该检验方法不被支持')
        if self.option("mann_ci") <= 0 or self.option("mann_ci") >= 1:
            raise OptionError('所输入的显著水平不在范围值内')
        if self.option("mann_type") not in ["double", "left", "right"]:
            raise OptionError('所输入的类型不在范围值内')
        if self.option("coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
            raise OptionError('秩和检验的置信区间的置信度不在范围值内')
        if self.option("compare_type") == "multi":
            if self.option("group_name").find(";") >= 0:
                g_list = self.option("group_name").split(';')
            elif self.option("group_name").find(",") >= 0:
                g_list = self.option("group_name").split(',')
            else:
                raise OptionError('分组名请按英文字符“;”或“,”分隔！')
            g_list = list(set(g_list))
            new_list = [i for i in g_list if i != '']
            if len(new_list) < 3:
                raise OptionError('至少选择三个不同分组进行比较')

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
        self.data_table_raw = self.work_dir + "/input_table_raw.txt"
        self.data_table = self.work_dir + "/input_table.txt"
        if not self.option("project_data"):
            if self.option("sample") != "column":
                with open(self.option("mann_input").prop["path"], 'r') as f:
                    sample_list = [i.rstrip().split('\t') for i in f.readlines()]
                    new_sample_list = map(lambda *a: '\t'.join(a) + '\n', *sample_list)
                    with open(self.data_table_raw, 'w') as fw1:
                        fw1.writelines(new_sample_list)
            else:
                if os.path.exists(self.data_table_raw):
                    os.remove(self.data_table_raw)
                os.link(self.option("mann_input").prop["path"], self.data_table_raw)
        self.name_dict = rename_name(self.data_table_raw,self.data_table)

    def run_mann(self):
        if self.option("compare_type") != "multi":
            if self.option("mann_type") == "double":
                mann_type = "two.side"
            elif self.option("mann_type") == "left":
                mann_type = "less"
            else:
                mann_type = "greater"
            options = {
                "mann_input": self.data_table,
                "mann_group": self.group_table,
                "sample": "column",
                "compare_type": self.option("compare_type"),
                "group_name1": self.option('group_name1'),
                "group_name2": self.option('group_name2'),
                "mann_ci": self.option('mann_ci'),
                "correction": self.option('correction'),
                "mann_type": mann_type,
                "coverage": self.option('coverage'),
                "method": self.option('method'),
            }
        else:
            options = {
                "mann_input": self.data_table,
                "mann_group": self.group_table,
                "sample": "column",
                "compare_type": self.option("compare_type"),
                "group_name": self.option('group_name'),
                "mann_ci": self.option('mann_ci'),
                "correction": self.option('correction'),
                "kru_methor": self.option('kru_methor'),
                "coverage": self.option('coverage'),
                "method": self.option('method'),
            }
        self.mann.set_options(options)
        self.mann.on("end", self.set_output, "mann")
        self.mann.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'mann':
            self.linkdir(obj.output_dir, 'mann')
        if os.path.exists(self.output_dir + "/mann/input_table.txt"):
            os.remove(self.output_dir + "/mann/input_table.txt")
        if os.path.exists(self.output_dir + "/mann/input_group.txt"):
            os.remove(self.output_dir + "/mann/input_group.txt")
        for file in os.listdir(self.output_dir + '/mann'):
            if file == "compare_result.xls":
                rename_name_back(self.output_dir + '/mann/compare_result.xls', self.name_dict, 1)
            else:
                rename_name_back(self.output_dir + '/mann/' + file, self.name_dict)
        if self.option("project_data"):
            if os.path.exists(self.output_dir +"/mann/input_data"):
                shutil.rmtree(self.output_dir +"/mann/input_data")
            os.mkdir(self.output_dir +"/mann/input_data")
            os.link(self.data_table_raw,self.output_dir +"/mann/input_data/input_table.txt")
            os.link(self.group_table, self.output_dir + "/mann/input_data/input_group.txt")
        for file in [self.output_dir+"/mann/compare_boxfile.xls",self.output_dir+"/mann/tmp_group_bar.xls"]:
            os.rename(file,file+"_tmp")
            with open(file+"_tmp") as f,open(file,"w") as t:
                title = f.readline()
                if not title.strip().startswith("Name"):
                    t.write("Name\t"+title)
                for x in f:
                    t.write(x)
            os.remove(file+"_tmp")
        self.set_db()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))



    def set_db(self):
        """
        保存结果到Mongo库
        """
        self.logger.info("保存结果到Mongo")
        main_id = self.option("main_id")
        mann_api = self.api.api("tool_lab.mann_test")
        mann_dir = self.output_dir + '/mann'
        box_file = self.output_dir + '/mann' + '/box_result_group.xls'
        if self.option("compare_type") != "multi":
            group_names = [self.option('group_name1'), self.option('group_name2')]
        else:
            if self.option("group_name").find(";") >= 0:
                group_names = self.option("group_name").split(';')
            else:
                group_names = self.option("group_name").split(',')
        mann_api.add_mann_detail(main_id, mann_dir, box_file, group_names,self.option('method'),self.option("compare_type"),self.option('kru_methor'),self.option('correction'))
        self.logger.info("保存结果到Mongo结束")
        self.end()

    def run(self):
        """
        运行
        """
        if self.option("project_data"):
            self.data_table_raw,self.group_table = meta_get_file(self.option('project_data'),self.work_dir)
        else:
            self.data_table_raw = self.option('mann_input').path
            self.group_table = self.option('mann_group').path
        self.file_rename()
        self.run_mann()
        super(MannTestWorkflow, self).run()

    def end(self):
        if self.option('mann_group'):
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "秩和检验分析结果输出目录"],
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "秩和检验分析结果输出目录"],
            ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MannTestWorkflow, self).end()
