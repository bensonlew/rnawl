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
import glob


class PlsdaRoplsWorkflow(Workflow):
    """

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PlsdaRoplsWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table", "type": "infile", "format": "tool_lab.table"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group"},  # 分组文件
            {"name": "groups","type": "string"},
            {"name": "group1","type": "string"},
            {"name": "group2","type": "string"},
            {'name': 'method', 'type': 'string', 'default': 'plsda'},  # 多元统计类型，pca，plsda, oplsda
            {'name': 'confidence', 'type': 'string', 'default': '0.95'},  # 置信度，与mul_type对应
            {'name': 'perm', 'type': 'string', 'default': '200'},  # 置换次数，与mul_type对应
            #{'name': 'data_trans', 'type': 'string', 'default': 'Par'},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.pls_tool = self.add_tool("tool_lab.diff_mul_stat")
        self.ellipse_tool = self.add_tool("tool_lab.ellipse")
        if self.option('method') == 'plsda':
            if '|' in self.option("groups"):
                self.group_names = self.option("groups").split('|')
            elif ',' in self.option("groups"):
                self.group_names = self.option("groups").split(",")
            elif '，' in self.option("groups"):
                self.group_names = self.option("groups").split("，")
            else:
                self.set_error('样本组名称必须用 ， 或 | 分割')
        else:
            self.group_names = [self.option("group1"), self.option("group2")]

    def sort_group_file_fun(self):
        group_file = self.option("group_table").prop["path"]
        group_map_sample = dict()
        with open(group_file) as f:
            f.readline()
            group_list = list()
            for line in f:
                line = line.strip()
                if line=='':
                    continue
                spline = line.split("\t")
                if spline[1] not in group_list:
                    group_list.append(spline[1])
                    group_map_sample[spline[1]] = list()
                group_map_sample[spline[1]].append(spline[0])
        group_file_sorted = self.work_dir + '/sorted_group.xls'
        with open(group_file_sorted, 'w') as fw:
            fw.write("#sample\tgroup_name\n")
            for g in group_list:
                #if g in  [self.option("group1"), self.option("group2")]: #self.option("groups").split(','):
                if g in self.group_names:
                    fw.write('\n'.join([i+'\t'+ g for i in sorted(group_map_sample[g])]) + '\n')
        return group_file_sorted

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        #if  self.option("group1") not in  self.option("group_table").prop['group_name_list'] or  self.option("group2") not in  self.option("group_table").prop['group_name_list']:
        for g in self.group_names:
            if g not in self.option("group_table").prop['group_name_list']:
                self.set_error('输入的组名没有在分组文件中')
        self.logger.info("start run  workflow")
        self.new_group_file = self.sort_group_file_fun()
        self.run_pls()
        super(PlsdaRoplsWorkflow, self).run()

    def run_pls(self):
        self.logger.info("module pls analysis start")
        options = {
            "exp_file": self.option("table"),
            "mul_type": self.option("method"),
            "group_file":  self.option("group_table").prop['path'],
            "groups" :   ','.join(self.group_names), #self.option("groups"),
            "confidence": self.option("confidence"),
            "perm": self.option("perm")
        }
        self.pls_tool.set_options(options)
        self.pls_tool.on('end', self.run_ellipse)
        self.pls_tool.run()

    def run_ellipse(self):
        self.logger.info("ellipse analysis start")
        pc_table = glob.glob(self.pls_tool.output_dir+'/*.sites.xls')[0]
        options = {
            "pc_table": pc_table,
            "group_table": self.new_group_file,
            "draw_all_sample" : 't'
        }
        self.ellipse_tool.set_options(options)
        self.ellipse_tool.on('end', self.set_db)
        self.ellipse_tool.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """

        pls_dir = self.pls_tool.output_dir
        self.move_file(pls_dir, self.output_dir)

        api_name = self.api.api("tool_lab.plsda_ropls")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！')

        api_name.add_exp_diff_comp(main_id, pls_dir, self.option("group_table").prop["path"])
        api_name.add_exp_diff_model(main_id, pls_dir)
        api_name.add_exp_diff_scatter(main_id, pls_dir)
        api_name.add_exp_diff_load(main_id,pls_dir)
        group_ci_ellipse = self.ellipse_tool.output_dir+'/ellipse_out.xls'
        all_sample_ellipse = self.ellipse_tool.work_dir+'/all_sample_ellipse.xls'
        api_name.add_ci_ellipse(main_id, group_ci_ellipse,all_sample_ellipse)

        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "样本PlsDA分析", 0]
        ]
        regexps = [
            [r".*.loadings\.xls", "xls", "主成分贡献度表", 0],
            [r".*.model\.xls", "xls", "模型参数表", 0],
            [r".*.permMN\.xls", "xls", "响应排序检验结果表", 0],
            [r".*.site\.xls", "xls", "样本各维度坐标", 0],
            [r".*.vip\.xls", "xls", "VIP值表", 0]
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(PlsdaRoplsWorkflow, self).end()


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
            "main_id" : "5e58bd0217b2bf4f565122e6",
            "group_table":"/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/remote_input/group_table/group.txt",
            "groups":"NG_D1,NG_D3,SMG_D1",
            "method":"plsda"
        }
    }

    wsheet = Sheet(data=data)

    wf = PlsdaRoplsWorkflow(wsheet)
    wf.run()