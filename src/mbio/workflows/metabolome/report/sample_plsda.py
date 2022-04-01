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


class SamplePlsdaWorkflow(Workflow):
    """
    多组 plsda
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SamplePlsdaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_abun"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_detail", "type": "string"},
            ###{'name': 'group_name', 'type': 'string', 'default': ''},  # 差异分组
            ###{'name': 'mul_type', 'type': 'string', 'default': 'pca;plsda;oplsda'},  # 多元统计类型，pca，plsda, oplsda
            {'name': 'confidence', 'type': 'string', 'default': '0.95'},  # 置信度，与mul_type对应
            {'name': 'perm', 'type': 'string', 'default': '200'},  # 置换次数，与mul_type对应
            {'name': 'data_trans', 'type': 'string', 'default': 'Par'},
            # 数据转化方法："UV","Ctr","Par"，""
            ###{'name': 'test_method', 'type': 'string', 'default': 't-test'},  # 差异检验方法
            ###{'name': 'side_type', 'type': 'string', 'default': 'two.side'},  # 单尾或双尾检验 two.side,less,greater
            {'name': 'table_type', 'type': 'string', 'default': 'pos'},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "test_api","type":"string","default":"f"}, #测试导表用
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.pls_tool = self.add_tool("metabolome.diff.diff_mul_stat")


    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run  workflow")
        if self.option("test_api") !='f':
            self.set_db()
        else:
            self.run_pls()
        super(SamplePlsdaWorkflow, self).run()

    def run_pls(self):
        """
        plsda
        """
        self.logger.info("module pls analysis start")
        options = {
            "exp_file": self.option("metab_table"),
            "metab_desc": self.option("metab_desc"),
            "mul_type": 'plsda',
            "group_file":  self.option("group"),
            ###"group_name": self.option("group_name"),
            "confidence": self.option("confidence"),
            "perm": self.option("perm"),
            "data_trans": self.option("data_trans"),
        }
        self.pls_tool.set_options(options)
        self.pls_tool.on('end', self.set_db)
        self.pls_tool.run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        if self.option("test_api") != 'f':
            pls_dir = self.option("test_api")
        else:
            pls_dir = self.pls_tool.output_dir
        self.move_file(pls_dir, self.output_dir)
        api_name = self.api.api("metabolome.sample_plsda")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！')

        api_name.add_sample_plsda(main_id,self.option("table_type"), main_id=main_id)
        api_name.add_exp_diff_bar(main_id, pls_dir)
        api_name.add_exp_diff_comp(main_id, pls_dir, self.option("group").prop["path"])
        api_name.add_exp_diff_model(main_id, pls_dir)
        api_name.add_exp_diff_scatter(main_id, pls_dir)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "expplsda",
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
            [".", "", "样本PlsDA分析", 0]
        ]
        regexps = [
            [r"PLS-DA\.loadings\.xls", "xls", "PLS-DA代谢物主成分贡献度表", 0],
            [r"PLS-DA\.model\.xls", "xls", "PLS-DA模型参数表", 0],
            [r"PLS-DA\.permMN\.xls", "xls", "PLS-DA响应排序检验结果表", 0],
            [r"PLS-DA\.sites\.xls", "xls", "PLS-DA样本各维度坐标", 0],
            [r"PLS-DA\.vips\.xls", "xls", "PLS-DA的VIP值表", 0]
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(SamplePlsdaWorkflow, self).end()


    def move_file(self, old_file, new_file):
        """
        递归移动文件夹的内容
        """
        if os.path.isfile(old_file):
            if not os.path.isdir(os.path.dirname(new_file)):
                os.makedirs(os.path.dirname(new_file))
            old_file_name = old_file.split("/")[-1]
            if not old_file_name  in [ "PLS-DA.ellipse.xls","PLS-DA.intercept.xls", "PLS-DA.loading.xls", "PLS-DA.vip.xls"]:
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
            self.set_error("链接失败：请检查%s", variables=(old_file), code="14701702")


if __name__ == '__main__':
    from biocluster.wsheet import Sheet

    data = {
        'name': 'test_sample_PLSDA',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "metab_table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            #"group_detail" : '{"a":["NG_D1_A","NG_D1_B"],"b":["NG_D1_C","NG_D1_D"]}',
            "main_table_id" : "5e58bd0217b2bf4f565122e5",
            "group":"/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/remote_input/group_table/group.txt",
            "metab_desc": "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_desc.txt",
            "test_api": "/mnt/ilustre/users/sanger-dev/workspace/20200311/SamplePlsda_tsg_36964/DiffMulStat/output"
        }
    }

    wsheet = Sheet(data=data)

    wf = SamplePlsdaWorkflow(wsheet)
    wf.run()