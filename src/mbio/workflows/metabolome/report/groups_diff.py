# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
import pandas as pd
import glob

class GroupsDiffWorkflow(Workflow):
    """
    多组差异分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GroupsDiffWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_abun"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "group_name", "type": "string"},
            #{'name': 'data_trans', 'type': 'string', 'default': 'UV;Par;Par'},
            # 数据转化方法："UV","Ctr","Par"，"", 与mul_type对应个数
            {"name": "test_method", "type" : "string", "default": "ow"},
            {'name': 'post_hoc', 'type': 'string', 'default': 'scheffe'},
            {'name': 'coverage', 'type': 'float', 'default':0.95},
            {'name': 'table_type', 'type': 'string', 'default': 'pos'}, ##'pos' or mix  or pos,neg
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "metab_table_neg","type": "infile", "format":"metabolome.metab_abun"},
            {"name": "metab_desc_neg", "type": "infile", "format": "sequence.profile_table"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_list = []
        self.tool_dic = {}

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")

        if self.option("metab_table_neg").is_set:
            self.run_groups_diff()
            self.run_groups_diff2()
            self.on_rely(self.tool_list,self.set_db)
            for i in self.tool_list:
                i.run()
        else:
            self.run_groups_diff()
            self.tool_list[0].on('end',self.set_db)
            self.tool_list[0].run()

        super(GroupsDiffWorkflow, self).run()


    def run_groups_fun(self,exp_file,desc_file):
        multiple = self.add_module("metabolome.groups_diff")
        options = {
            "metab_table": exp_file,
            "group" :  self.option("group"),
            "metab_desc" : desc_file,
            "group_name" : self.option("group_name"),
            "test_method" : self.option("test_method"),
            "post_hoc"  : self.option("post_hoc"),
            "coverage" : self.option("coverage")
        }
        multiple.set_options(options)
        return multiple

    def run_groups_diff(self):
        self.logger.info("start run group diff!")
        diff_module = self.run_groups_fun(self.option("metab_table").path, self.option("metab_desc").path)
        self.tool_list.append(diff_module)
        self.tool_dic['not_neg'] = diff_module

    def run_groups_diff2(self):
        self.logger.info("start run group diff 2 !")
        diff_module = self.run_groups_fun(self.option("metab_table_neg").path, self.option("metab_desc_neg").path)
        self.tool_list.append(diff_module)
        self.tool_dic['neg'] = diff_module

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.groups_diff")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14701701")

        if len(self.tool_list) == 1:
            tool = self.tool_list[0]
            type = self.option("table_type")  #'pos' or mix
            #link to output
            target_dir = self.output_dir + '/' + type
            if os.path.exists(target_dir):
                shutil.rmtree(target_dir)
            shutil.copytree(tool.output_dir, target_dir)
            #api
            api_name.add_groups_diff_detail(target_dir,main_id,table_type=type)

        elif len(self.tool_list) == 2:
            for k in self.tool_dic.keys():
                tool = self.tool_dic[k]
                if k == 'not_neg':
                    type = 'pos'
                else:
                    type = k
                #link to output
                target_dir = self.output_dir + '/' + type
                if os.path.exists(target_dir):
                    shutil.rmtree(target_dir)
                shutil.copytree(tool.output_dir, target_dir)
                #api
                api_name.add_groups_diff_detail(target_dir,main_id,table_type=type)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "diff_multigroup",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "多组比较统计结果文件夹", 0 ],
            ["./pos", "", "差异检验结果", 0],
            ["./neg", "", "差异检验结果", 0],
            ["./mix", "", "差异检验结果", 0],
        ]
        regexps = [
            [r".*/.*.xls", "xls", "结果表", 0]
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(GroupsDiffWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test_sample_diff',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "metab_table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "metab_desc" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_desc.txt",
            "group_detail" : '{"a":["NG_D1_A"],"b":["NG_D1_D"]}',
            "main_table_id" : "5e4ce5e817b2bf4b326f4b20",
            "group_name" : "test",
            "group" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/remote_input/group_table/group.txt",
            "metab_table_neg" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "metab_desc_neg" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_desc.txt"
        }
    }

    wsheet = Sheet(data=data)

    wf = GroupsDiffWorkflow(wsheet)
    wf.run()