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


class PlsdaMixomicsWorkflow(Workflow):
    """

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PlsdaMixomicsWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table", "type": "infile", "format": "tool_lab.table"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group"},  # 分组文件
            {"name": "groups","type": "string"},
            #{"name": "group1","type": "string"},
            #{"name": "group2","type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

        self.pls_tool = self.add_tool("tool_lab.plsda")
        self.ellipse = self.add_tool("tool_lab.ellipse")
        if '|' in self.option("groups"):
            self.group_names = self.option("groups").split('|')
        elif ',' in self.option("groups"):
            self.group_names = self.option("groups").split(",")
        elif '，' in self.option("groups"):
            self.group_names = self.option("groups").split("，")
        else:
            self.set_error('样本组名称必须用 ， 或 | 分割')

    def pick_data(self):
        group_data = pd.read_table(self.option("group_table").path,sep='\t')
        c0 = group_data.columns[0]
        c1 = group_data.columns[1]
        #pick_group = self.option("groups").split(',')
        pick_group = self.group_names
        #pick_group = [self.option('group1'), self.option('group2')]
        group_data = group_data[group_data[c1].map(lambda x: x in pick_group)]
        self.new_group = self.option("group_table").path + '_new'
        group_data.to_csv(self.new_group, sep='\t', index=False)

        #sample_list = group_data[c0].tolist()
        #table_data = pd.read_table(self.option("table").path, sep='\t')
        self.table = self.option("table").path




    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        for group_name in self.group_names:
            if group_name not in self.option("group_table").prop['group_name_list']:
                self.set_error("输入的组名没有在分组文件中")

        self.logger.info("start run  workflow")
        self.run_pls()
        super(PlsdaMixomicsWorkflow, self).run()

    def run_pls(self):
        self.pick_data()
        self.logger.info("module pls analysis start")
        options = {
            "otutable": self.table,
            "group":  self.new_group,
            "ellipse" : 'T'
        }
        self.pls_tool.set_options(options)
        self.pls_tool.on('end', self.run_ellipse)
        self.pls_tool.run()


    def run_ellipse(self):
        options = {}
        options['group_table'] = self.new_group
        options['analysis'] = 'plsda'
        options['meta'] = ''
        options['pc_table'] = self.pls_tool.output_dir + '/plsda_sites.xls'
        options['draw_all_sample'] = 't'
        self.ellipse.set_options(options)
        self.ellipse.on('end', self.set_db)
        self.ellipse.run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """

        pls_dir = self.pls_tool.output_dir
        self.move_file(pls_dir, self.output_dir)
        api_name = self.api.api("tool_lab.plsda_mixomics")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        api_name.add_plsda_result(pls_dir, main_id, self.option("group_table").prop["path"])
        api_name.insert_plsda_meta_comp_ellipse(self.ellipse.output_dir+'/ellipse_out.xls', main_id,self.ellipse.work_dir+'/all_sample_ellipse.xls')
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "样本PlsDA分析", 0],
            ["plsda_sites.xls", "xls", "样本坐标表", 0],
            ["plsda_rotation.xls", "xls", "物种主成分贡献度表", 0],
            ["plsda_importance.xls", "xls", "主成分组别特征值表", 0],
            ["plsda_importancepre.xls", "xls", "主成分解释度表", 0]
        ]

        result_dir.add_relpath_rules(relpath_rules)
        super(PlsdaMixomicsWorkflow, self).end()


    def move_file(self, old_file, new_file):
        """
        递归移动文件夹的内容
        """
        if os.path.isfile(old_file):
            if not os.path.isdir(os.path.dirname(new_file)):
                os.makedirs(os.path.dirname(new_file))
            if os.path.exists(new_file):
                os.remove(new_file)
            os.link(old_file, new_file)
        elif os.path.isdir(old_file):
            if not os.path.exists(new_file):
                os.makedirs(new_file)
            for file in os.listdir(old_file):
                file_path = os.path.join(old_file, file)
                new_path = os.path.join(new_file, file)
                self.move_file(file_path, new_path)
        else:
            self.set_error("链接失败：请检查%s", variables=(old_file))


if __name__ == '__main__':
    from biocluster.wsheet import Sheet

    data = {
        'name': 'test_sample_PLSDA',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "main_id" : "5e58bd0217b2bf4f565123e9",
            "group_table":"/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/remote_input/group_table/group.txt",
            "groups":"NG_D1,NG_D3,SMG_D1"
            #"group1" : ''
        }
    }

    wsheet = Sheet(data=data)

    wf = PlsdaMixomicsWorkflow(wsheet)
    wf.run()