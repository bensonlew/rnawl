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


class ExpDiffWorkflow(Workflow):
    """
    代谢样本相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpDiffWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_abun"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_detail", "type": "string"},
            {'name': 'group_name', 'type': 'string', 'default': ''},  # 差异分组
            {'name': 'mul_type', 'type': 'string', 'default': 'pca;plsda;oplsda'},  # 多元统计类型，pca，plsda, oplsda
            {'name': 'confidence', 'type': 'string', 'default': '0.95;0.95;0.95'},  # 置信度，与mul_type对应
            {'name': 'perm', 'type': 'string', 'default': '0;200;200'},  # 置换次数，与mul_type对应
            {'name': 'data_trans', 'type': 'string', 'default': 'UV;Par;Par'},
            # 数据转化方法："UV","Ctr","Par"，"", 与mul_type对应个数
            {'name': 'test_method', 'type': 'string', 'default': 't-test'},  # 差异检验方法
            {'name': 'side_type', 'type': 'string', 'default': 'two.side'},  # 单尾或双尾检验 two.side,less,greater
            {'name': 'table_type', 'type': 'string', 'default': 'pos'},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.diff_pls = self.add_module("metabolome.diff_pls")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Exp_diff workflow")
        self.run_diff_pls()
        super(ExpDiffWorkflow, self).run()

    def run_diff_pls(self):
        self.logger.info("start run diff_pls !")
        group_name = self.option("group_name")
        group_table = self.option("group")
        self.logger.info("------------------")
        self.logger.info(self.option("perm"))
        self.logger.info(self.option("test_method"))
        self.logger.info(self.option("side_type"))
        exp_file = self.option("metab_table")
        options = {
            'exp_file': exp_file,
            "group_file": group_table,
            'group_name': group_name,
            'mul_type': 'pca;plsda;oplsda',
            'confidence': self.option("confidence"),
            'perm': self.option("perm"),
            'data_trans': self.option("data_trans"),
            'test_method': self.option("test_method"),
            'side_type': self.option("side_type"),
            'metab_desc': self.option("metab_desc")
        }
        self.diff_pls.set_options(options)
        self.diff_pls.on('end', self.set_db)
        self.diff_pls.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.exp_diff")
        self.link_pls()
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14701701")
        sanger_type, sanger_path = self._sheet.output.split(':')
        my_bucket = Config().get_project_region_bucket(project_type="metabolome")
        if sanger_type in my_bucket:
            #web_path =  sanger_path.lstrip("//") + "DiffTest/"
            web_path = os.path.join(self._sheet.output,"DiffTest/")
        else:
            web_path =  sanger_path + "/DiffTest/"
        #diff_file_dir = self.diff_pls.option("diffStat_dir").prop["path"]
        id_diff_file_dir = self.diff_pls.option("id_diffStat_dir").prop["path"]
        pls_dir = self.diff_pls.option("pls_dir").prop["path"]
        api_name.add_exp_diff(main_id,self.option("table_type"), main_id=main_id, diff_dir=web_path)
        api_name.add_exp_diff_detail(main_id, id_diff_file_dir, self.option("table_type"))
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
                "submit_loc": "diff_twogroup",
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
            [".", "", "差异代谢物计算与统计结果文件夹", 0, "150048"],
            ["./DiffMulStat", "", "多元统计模型结果", 0, "150050"],
            ["./DiffTest", "", "差异检验结果", 0, "150064"]
        ]
        regexps = [
            [r"DiffMulStat/.*_vs_.*/OPLS-DA\.model\.xls", "xls", "OPLS-DA模型参数表", 0, "150051"],
            [r"DiffMulStat/.*_vs_.*/OPLS-DA\.permMN\.xls", "xls", "OPLS-DA响应排序检验结果表", 0, "150052"],
            [r"DiffMulStat/.*_vs_.*/OPLS-DA\.loadings\.xls", "xls", "OPLS-DA代谢物主成分贡献度表", 0, "150053"],
            [r"DiffMulStat/.*_vs_.*/OPLS-DA\.sites\.xls", "xls", "OPLS-DA样本各维度坐标", 0, "150054"],
            [r"DiffMulStat/.*_vs_.*/OPLS-DA\.vips\.xls", "xls", "OPLS-DA的VIP值表", 0, "150055"],
            [r"DiffMulStat/.*_vs_.*/PCA\.loadings\.xls", "xls", "PCA代谢物主成分贡献度表", 0, "150014"],
            [r"DiffMulStat/.*_vs_.*/PCA\.model\.xls", "xls", "PCA模型参数表", 0, "150015"],
            [r"DiffMulStat/.*_vs_.*/PCA\.sites\.xls", "xls", "PCA样本各维度坐标", 0, "150016"],
            [r"DiffMulStat/.*_vs_.*/PLS-DA\.loadings\.xls", "xls", "PLS-DA代谢物主成分贡献度表", 0, "150056"],
            [r"DiffMulStat/.*_vs_.*/PLS-DA\.model\.xls", "xls", "PLS-DA模型参数表", 0, "150057"],
            [r"DiffMulStat/.*_vs_.*/PLS-DA\.permMN\.xls", "xls", "PLS-DA响应排序检验结果表", 0, "150058"],
            [r"DiffMulStat/.*_vs_.*/PLS-DA\.sites\.xls", "xls", "PLS-DA样本各维度坐标", 0, "150059"],
            [r"DiffMulStat/.*_vs_.*/PLS-DA\.vips\.xls", "xls", "PLS-DA的VIP值表", 0, "150060"],
            [r"DiffTest/.*_vs_.*\.diff\.exp\.xls", "xls", "差异检验结果表", 0, "150061"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(ExpDiffWorkflow, self).end()

    def link_pls(self):
        all_dirs = os.listdir(self.diff_pls.work_dir)
        newdir = self.output_dir
        oldfiles = []
        newfiles = []
        for eachdir in all_dirs:
            if eachdir == "DiffMulStat":
                eachdirpath = self.diff_pls.work_dir + "/DiffMulStat/output"
                allfiles = os.listdir(eachdirpath)
                newdir = os.path.join(self.output_dir, eachdir)
                for i in allfiles:
                    oldfile = os.path.join(eachdirpath, i)
                    oldfiles.append(oldfile)
                    newfile = os.path.join(newdir, i)
                    newfiles.append(newfile)
            if eachdir == "MergeDiff":
                DiffTestdir = self.diff_pls.work_dir + "/MergeDiff/output/DiffStat"
                allfiles = os.listdir(DiffTestdir)
                newdir = newdir = os.path.join(self.output_dir, "DiffTest")
                for i in allfiles:
                    oldfile = os.path.join(DiffTestdir, i)
                    oldfiles.append(oldfile)
                    newfile = os.path.join(newdir, i)
                    newfiles.append(newfile)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(oldfiles)):
            self.move_file(oldfiles[i], newfiles[i])


    def move_file(self, old_file, new_file):
        """
        递归移动文件夹的内容
        """
        if os.path.isfile(old_file):
            if not os.path.isdir(os.path.dirname(new_file)):
                os.makedirs(os.path.dirname(new_file))
            old_file_name = old_file.split("/")[-1]
            if not old_file_name  in ["PCA.ellipse.xls", "OPLS-DA.ellipse.xls", "OPLS-DA.intercept.xls",
                         "OPLS-DA.loading.xls", "OPLS-DA.vip.xls", "PCA.loading.xls", "PLS-DA.ellipse.xls",
                         "PLS-DA.intercept.xls", "PLS-DA.loading.xls", "PLS-DA.vip.xls"]:
                os.link(old_file, new_file)
        elif os.path.isdir(old_file):
            os.makedirs(new_file)
            for file in os.listdir(old_file):
                file_path = os.path.join(old_file, file)
                new_path = os.path.join(new_file, file)
                self.move_file(file_path, new_path)
        else:
            self.set_error("链接失败：请检查%s", variables=(old_file), code="14701702")
