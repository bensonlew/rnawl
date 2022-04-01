# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0525

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types



class TwoSampleWorkflow(Workflow):
    """
    代谢样本相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TwoSampleWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_abun"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_detail", "type": "string"},
            {'name': 'diff_group_name', 'type': 'string', 'default': ''},  # 差异分组
            #{'name': 'data_trans', 'type': 'string', 'default': 'UV;Par;Par'},
            # 数据转化方法："UV","Ctr","Par"，"", 与mul_type对应个数
            {'name': 'test_method', 'type': 'string', 'default': 'fisher'},  # 差异检验方法  "chiq", 'fisher'
            {'name': 'tail', 'type': 'string', 'default': 'two-tailed'},  # # "two-tailed", "left-tailed", "right-tailed"
            {'name': 'correct', 'type': 'string', 'default':'bonferroni'},
            {'name': 'table_type', 'type': 'string', 'default': 'pos'}, ##'pos' or mix  or pos,neg
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "metab_table_neg","type": "infile", "format":"metabolome.metab_abun"},
            {"name": "metab_desc_neg", "type": "infile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_list = []
        self.tool_dic = {}

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")
        tail_change = {
            'two-tailed':'two.side',
            'left-tailed': 'less',
            'right-tailed': 'greater'
        }

        self.tail = tail_change[self.option('tail')]

        if self.option("metab_table_neg").is_set:
            self.run_diff_sample()
            self.run_diff_sample2()
            self.on_rely(self.tool_list,self.set_db)
            for i in self.tool_list:
                i.run()
        else:
            self.run_diff_sample()
            self.tool_list[0].on('end',self.set_db)
            self.tool_list[0].run()

        super(TwoSampleWorkflow, self).run()


    def run_diff_sample(self):
        self.logger.info("start run diff_sample!")
        diff_module = self.add_module('metabolome.two_sample_batch')

        options = {
            "metab_table": self.option('metab_table'),
            "group_detail" : self.option("group_detail"),
            "diff_group_name" : self.option('diff_group_name'),
            "test_method" : self.option("test_method"),
            "side_type" : self.tail,
            "correct" : self.option("correct"),
            'metab_desc': self.option("metab_desc")
        }
        diff_module.set_options(options)
        self.tool_list.append(diff_module)
        self.tool_dic['not_neg'] = diff_module

    def run_diff_sample2(self):
        self.logger.info("start run diff_sample 2 !")
        diff_module = self.add_module('metabolome.two_sample_batch')
        options = {
            "metab_table": self.option("metab_table_neg"),
            "group_detail" : self.option("group_detail"),
            "diff_group_name" : self.option('diff_group_name'),
            "test_method" : self.option("test_method"),
            "side_type" : self.tail,
            "correct" : self.option("correct"),
            'metab_desc': self.option("metab_desc_neg")
        }
        diff_module.set_options(options)
        self.tool_list.append(diff_module)
        self.tool_dic['neg'] = diff_module


    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.two_sample")
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
            api_name.add_two_sample_detail(tool.output_dir,main_id,table_type=type)
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
                api_name.add_two_sample_detail(tool.output_dir,main_id,table_type=type)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "差异代谢物计算与统计结果文件夹", 0, ],
            ["./pos", "", "差异检验结果", 0],
            ["./neg", "", "差异检验结果", 0],
            ["./mix", "", "差异检验结果", 0],
        ]
        regexps = [
            [r".*/.*_vs_.*/.*_result.xls", "xls", "p值表", 0],
            [r".*/.*_vs_.*/.*_CI.xls", "xls", "CI表", 0]
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(TwoSampleWorkflow, self).end()


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
            "diff_group_name" : "a_vs_b;b_vs_a" ,
            "metab_table_neg" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "metab_desc_neg" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_desc.txt"
        }
    }

    wsheet = Sheet(data=data)

    wf = TwoSampleWorkflow(wsheet)
    wf.run()