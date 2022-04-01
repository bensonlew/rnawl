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
            {"name": "protein_table", "type": "infile", "format": "metabolome.metab_abun"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            # {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_dict", "type": "string"},
            {"name": "group_name", "type": "string"},
            #{'name': 'data_trans', 'type': 'string', 'default': 'UV;Par;Par'},
            # 数据转化方法："UV","Ctr","Par"，"", 与mul_type对应个数
            {"name": "test_method", "type" : "string", "default": "ow"},
            {'name': 'post_hoc', 'type': 'string', 'default': 'scheffe'},
            {'name': 'coverage', 'type': 'float', 'default':0.95},
            # {'name': 'table_type', 'type': 'string', 'default': 'pos'}, ##'pos' or mix  or pos,neg
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            # {"name": "protein_table_neg","type": "infile", "format":"metabolome.metab_abun"},
            # {"name": "metab_desc_neg", "type": "infile", "format": "sequence.profile_table"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_list = []

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")

        self.run_groups_diff()
        self.tool_list[0].on('end',self.set_db)
        self.tool_list[0].run()

        super(GroupsDiffWorkflow, self).run()


    def run_groups_fun(self,exp_file):
        multiple = self.add_module("itraq_and_tmt.groups_diff")
        options = {
            "protein_table": exp_file,
            "group" :  self.option("group"),
            # "metab_desc" : desc_file,
            "group_name" : self.option("group_name"),
            "test_method" : self.option("test_method"),
            "post_hoc"  : self.option("post_hoc"),
            "coverage" : self.option("coverage")
        }
        multiple.set_options(options)
        return multiple

    def run_groups_diff(self):
        self.logger.info("start run group diff!")
        diff_module = self.run_groups_fun(self.option("protein_table").path)
        self.tool_list.append(diff_module)

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("itraq_and_tmt.groups_diff")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14701701")

        tool = self.tool_list[0]
        # type = self.option("table_type")  #'pos' or mix
        #link to output
        # shutil.copytree(tool.output_dir, self.output_dir)
        #api
        api_name.add_groups_diff_detail(tool.output_dir,main_id)
        self.end()

    def end(self):
        tool = self.tool_list[0]
        result_dir = self.add_upload_dir(tool.output_dir)
        relpath_rules = [
            [".", "", "多组比较统计结果文件夹", 0 ],
        ]
        regexps = [
            [r".*.xls", "xls", "结果表", 0]
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
            "protein_table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "metab_desc" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_desc.txt",
            "group_dict" : '{"a":["NG_D1_A"],"b":["NG_D1_D"]}',
            "main_table_id" : "5e4ce5e817b2bf4b326f4b20",
            "group_name" : "test",
            "group" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/remote_input/group_table/group.txt",
            "protein_table_neg" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "metab_desc_neg" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_desc.txt"
        }
    }

    wsheet = Sheet(data=data)

    wf = GroupsDiffWorkflow(wsheet)
    wf.run()