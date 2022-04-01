# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0709

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


class MetabsetVipWorkflow(Workflow):
    """
    代谢样本相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetVipWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'diff_dir', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},  # 差异分析结果目录
            {'name': 'diff_dir2', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "group_detail", "type": "string"},
            {"name": "table_type", "type": "string"},
            #{'name': 'group_name', 'type': 'string', 'default': ''},  # 使用的分组ing，“P_vs_C;C_vs_X”
            {'name': 'group_method', 'type': 'string', 'default': 'none'},  # 分组计算方法
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_set_id", "type": "string"},
            {'name': 'vip_type', 'type': 'string', 'default': 'oplsda'},  # vip来源 plsda，oplsda
            {'name': 'vip_cut', 'type': 'float', 'default': 1.0},  # vip 阈值
            {'name': 'vip_top', 'type': 'int', 'default': 30},  # Top vip 数目
            {'name': 'mct', 'type': 'string', 'default': ''},  # 代谢物聚类算法，hierarchy、kmeans、无
            {'name': 'mcd', 'type': 'string', 'default': ''},  # 代谢物距离计算方式
            {'name': 'n_cluster', 'type': 'int', 'default': 0},  # 代谢物聚类数目，kmeans时使用
            {'name': 'mcm', 'type': 'string', 'default': ''},  # 代谢物聚类方式, hierarchy时使用，"complete","average","single"
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "scale", "type": "string", "default": "scale"},  ## 是否标准化聚类数据，scale, unscale
            {"name": "metab_table", "type": "infile", "format": "sequence.profile_table"}, # 原始所有样本丰度表，scale时使用
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.metab_vip = self.add_tool("metabolome.metabset.metab_vip")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Vip workflow")
        self.run_vip()
        super(MetabsetVipWorkflow, self).run()

    def run_vip(self):
        self.logger.info("start run_vip !")
        diff_dir = self.option("diff_dir").prop["path"]

        # 20190806 v2
        if self.option('diff_dir2').is_set:
            diff_dir2 = self.option('diff_dir2').prop['path']
            all_dir = self.work_dir + '/all_diff_dir'
            if not os.path.exists(all_dir):
                os.mkdir(all_dir)
            files = os.listdir(diff_dir)
            files2 = os.listdir(diff_dir2)
            #base_dir = os.path.dirname(diff_dir)
            #base_dir2 = os.path.dirname(diff_dir2)
            share_file = set(files) & set(files2)
            for f in share_file:
                d1 = pd.read_table(diff_dir+'/'+f,header=0)
                d2 = pd.read_table(diff_dir2+'/'+f,header=0)
                data_cat = pd.concat([d1,d2],axis=0)  #LC 阴阳离子表合并
                data_cat.to_csv(all_dir+'/'+f,sep='\t' ,index=False)
            diff_dir = all_dir

        options = {
            'diff_dir': diff_dir,
            'metab_trans': self.option("metab_desc"),
            'vip_type': self.option("vip_type"),
            'vip_cut': self.option("vip_cut"),
            'vip_top': self.option("vip_top"),
            'mct': self.option("mct")
        }
        if self.option("mct") == "hierarchy":
            options["mcd"] = self.option("mcd")
            options["mcm"] = self.option("mcm")
        if self.option("mct") == "kmeans":
            options["mcd"] = self.option("mcd")
            options["n_cluster"] = self.option("n_cluster")
        if self.option("group_method") != "none":
            options["group_method"] = self.option("group_method")
            options["group"] = self.option("group")
        if self.option("metab_set_table").is_set:
            options["metab_set_table"] = self.option("metab_set_table")
        if "scale" in self.get_option_object().keys() and self.option("scale") == "scale":
            options["scale"] = True
            options["metab_table"] = self.option("metab_table")
        #if self.option("group_name"):
        #    options["group_name"] = self.option("group_name")
        self.metab_vip.set_options(options)
        self.metab_vip.on('end', self.set_db)
        self.metab_vip.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.metabset_vip")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14701801")
        vip_dir = self.metab_vip.output_dir
        self.logger.info(vip_dir)
        api_name.add_metabset_vip(main_id, main_id=main_id)
        #table_type = self.option("table_type")
        api_name.add_metabset_vip_detail(main_id, vip_dir, self.option("vip_type"), self.option('metab_desc').path,
                                         scale=self.option("scale"),out_dir=vip_dir)
        self.output_dir = self.metab_vip.output_dir
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetvip",
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
            [".", "", "VIP分析结果文件夹", 0, "150040"]
        ]
        regexps = [
            [r".*/metab.kmeans_cluster.xls", "", "代谢物kmeans分类文件", 0, "150065"],
            [r".*/metab.cluster_tree.xls", "xls", "代谢物树文件", 0, "150028"],
            [r".*/Vip_exp.xls", "xls", "VIP值表", 0, "150049"],
            [r".*/metab.subcluster_.*\.xls", "xls", "代谢物kmeans各子类结果表", 0, "150066"],
            [r".*/metab_id.cluster_tree.xls", "xls", "代谢物id聚类树文件", 0, "150075"],
            [r".*/Vip_scale_exp.xls", "xls", "代谢物VIP值与标准化后表达量", 0, "150078"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(MetabsetVipWorkflow, self).end()
