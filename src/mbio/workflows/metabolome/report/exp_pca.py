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


class ExpPcaWorkflow(Workflow):
    """
    代谢pca分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpPcaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "group_detail", "type": "string"},
            {"name": 'mul_type', "type": "string", "default": "pca"},
            {'name': 'confidence', 'type': 'string', 'default': '0.95'},
            {"name": "data_trans", "type": "string", "default": "UV"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.pca = self.add_tool("metabolome.diff.diff_mul_stat")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Exp_pca workflow")
        self.run_pca()
        super(ExpPcaWorkflow, self).run()

    def run_pca(self):
        self.logger.info("start run pca !")
        profile = self.option("metab_table").prop["path"]
        self.group = self.option("group_table").prop["path"]
        options = {
            'exp_file': profile,
            "group_file": self.group,
            'mul_type': self.option("mul_type"),
            'confidence': self.option("confidence"),
            'data_trans': self.option("data_trans"),
            'metab_desc': self.option("metab_desc")
        }
        self.pca.set_options(options)
        self.pca.on('end', self.set_db)
        self.pca.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.exp_pca")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14701401")
        loadings = self.link_file("PCA.loadings.xls", "PCA.loadings.xls")
        model = self.link_file("PCA.model.xls", "PCA.model.xls")
        sites = self.link_file("PCA.sites.xls", "PCA.sites.xls")
        ellipse = os.path.join(self.pca.output_dir, "PCA.ellipse.xls")
        api_name.add_exp_pca(main_id=main_id)
        api_name.add_exp_pca_detail(main_id, sites, ellipse, group_file=self.group)
        api_name.add_exp_pca_model(main_id, model)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "exppca",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.pca.output_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)
        return newfile

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "PCA结果文件夹", 0, "150013"],
            ["PCA.loadings.xls", "xls", "PCA代谢物主成分贡献度表", 0, "150014"],
            ["PCA.model.xls", "xls", "PCA模型参数表", 0, "150015"],
            ["PCA.sites.xls", "xls", "PCA样本各维度坐标", 0, "150016"],
        ])
        super(ExpPcaWorkflow, self).end()
