# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import os,shutil
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
from biocluster.core.exceptions import OptionError
import datetime
from biocluster.file import download
from collections import defaultdict
import shutil
import gevent
from mbio.packages.bac_comp_genome.common_function import link_dir


class BetaDiversityWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        细菌比较基因组泛基因组分析的第一步：
        聚类、泛基因组大小计算和变异分析；
        第二步 对聚类结果进行分组方案的合并
        第三步 对聚类的结果进行venn图的分析
        :return:
        """
        self._sheet = wsheet_object
        super(BetaDiversityWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "analysis", "type": "string","default": "pca"},  ### 分析名称
            {"name": "abundance_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"}, ### 丰度表
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"}, ###group表
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"}, ###距离矩阵
            {"name": "dis_method", "type": "string", "default": ""}, ###计算距离矩阵的方法
            {"name": "hcluster_method", "type": "string", "default": ""}, ###层级聚类的方法
            {"name": "sample_method", "type": "string", "default": ""},  # 样本的合并方式, ""为不进行合并
            {'name': 'corr_method', 'type': 'string', 'default': 'pearsonr'},
            {'name': 'function', 'type': 'string'}, ##对应功能水平的名称
            {'name': 'group_detail', 'type': 'string'}, ##to_file专用
            {'name': 'file_path', 'type': 'string'} , ##对应cog和kegg注释结果计算的丰度表
            {"name": "ellipse", "type": "string", "default": "T"},
            {"name": "is_group", "type": "string", "default": "no"},
            {"name": "group_method", "type": "string", "default": ""} ## 这里的group_method有四种["", "sum", "average", "middle"]
        ]
        self.beta_diversity = self.add_module("bac_comp_genome.beta_diversity")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.function = {
            "cog" : "all.cog_abundance.xls",
            "function" : "all.Function_abundance.xls",
            "ko" : "all.KO_abundance.xls",
            "pathway" : "all.level3_abundance.xls",
        }

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        self.logger.info("开始pan的参数检查")
        if not self.option("analysis"):
            raise OptionError("必须设置输入的analysis")
        if not (self.option("function") and self.option("corr_method")):
            raise OptionError("请提供输入丰度文件！")
        if not self.option("group").is_set:
            raise OptionError("请提供输入group文件！")
        if not self.option("main_id"):
            raise OptionError("请提供输入main_id")
        return True

    def download_files(self):
        """
        根据功能下载得到file文件
        :return:
        """
        file_path = self.option("file_path")
        file_name = self.function[self.option("function")]
        download_file_path = os.path.join(file_path, file_name)
        self.file_path = os.path.join(self.work_dir, 'abundance_table.xls')
        if download_file_path.startswith("/mnt/"):
            os.system("cp {} {}".format(download_file_path, self.file_path))
        else:
            if os.path.exists(self.file_path):
                os.remove(self.file_path)
            download(download_file_path, self.file_path)

    def run_beta_diversity(self):
        """
        基因组聚类分析
        :return:
        """
        self.logger.info("开始进行差异显著性分析")

        opts = ({
            "analysis": self.option("analysis"),
            "abundance_table": self.file_path
            })
        # if self.option("abundance_table").is_set:
        #     opts["abundance_table"] = self.option("abundance_table")
        # else:
        #     opts['abundance_table'] = self.file_path
        if self.option("group").is_set:
            opts["group"] = self.option("group")
        if self.option("dis_method") != '':
            opts["dis_method"] = self.option("dis_method")
        if self.option("hcluster_method") != '':
            opts["hcluster_method"] = self.option("hcluster_method")
        if self.option("sample_method") != '':
            opts["sample_method"] = self.option("sample_method")
        if self.option("corr_method") != '':
            opts["corr_method"] = self.option("corr_method")
        if self.option("ellipse") in ['T', 'F']:
            opts["ellipse"] = self.option("ellipse")
        if self.option("is_group") in ['yes', 'no']:
            opts["is_group"] = self.option("is_group")
        if self.option("group_method") in ['', 'sum', "average", "middle"]:
            opts["group_method"] = self.option("group_method")
        self.beta_diversity.set_options(opts)
        self.beta_diversity.run()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('设置结果目录')
        result_dir = self.beta_diversity.output_dir
        link_dir(result_dir, self.output_dir)
        self.logger.info('设置结果目录成功')
        self.set_db()

    def set_db(self):
        """
        将数据导入mongo
        :return:
        """
        self.logger.info('正在写入mongo数据库')
        # self.remote_dir = self._sheet.output
        self.remote_dir = self.output_dir
        main_id = self.option("main_id")
        if self.option("analysis") in ["pca"]:
            web_path = os.path.join(self.remote_dir, "Pca")
            data_path = os.path.join(self.output_dir, "Pca")
            api_pca_path = self.api.api("bac_comp_genome.pca")
            api_pca_path.add_pca_data(main_id, web_path, pca_data=data_path)
        elif self.option("analysis") in ["pcoa"]:
            api_pcoa_path = self.api.api("bac_comp_genome.pcoa")
            data_path = os.path.join(self.output_dir, "Pcoa")
            api_pcoa_path.add_pcoa_data(main_id, pcoa_data=data_path)
        elif self.option("analysis") in ["nmds"]:
            api_nmds_path = self.api.api("bac_comp_genome.nmds")
            data_path = os.path.join(self.output_dir, "Nmds")
            api_nmds_path.add_nmds_data(main_id, nmds_data=data_path)
        elif self.option("analysis") in ["hcluster"]:
            api_hcluster_path = self.api.api("bac_comp_genome.hcluster")
            data_path = os.path.join(self.output_dir, "Hcluster")
            api_hcluster_path.add_hcluster_data(main_id, hcluster_data=data_path)
        elif self.option("analysis") in ["correlation"]:
            api_correlation_path = self.api.api("bac_comp_genome.correlation")
            data_path = os.path.join(self.output_dir, "Correlation")
            api_correlation_path.add_correlation_data(main_id, correlation_data=data_path, method=self.option("corr_method"))
        else:
            self.set_error("不存在此分析内容")
        self.end()

    def run(self):
        """
        开始运行了
        :return:
        """
        self.logger.info("开始运行")
        self.beta_diversity.on('end', self.set_output)
        self.download_files()
        self.run_beta_diversity()
        super(BetaDiversityWorkflow, self).run()

    def end(self):
        """
        结束了
        :return:
        """
        repaths = [
            ["Pca", "", "PCA分析结果输出目录",0,""],
            ["Pca/pca_sites.xls", "xls", "PCA样本各成分轴坐标",0,""],
            ["Pca/pca_rotation.xls", "xls", "PCA主成分贡献度表",0,""],
            ["Pca/pca_importance.xls", "xls", "PCA主成分解释度表",0,""],
            ["Pca/pca_rotation_all.xls", "xls", "PCA全部主成分贡献度表",0,""],
            ["Pcoa", "", "PCoA分析结果目录", 0, ""],
            ["Pcoa/pcoa_eigenvalues.xls", "xls", "PCoA矩阵特征值", 0, ""],
            ["Pcoa/pcoa_eigenvaluespre.xls", "xls", "PCoA特征解释度百分比", 0, ""],
            ["Pcoa/pcoa_sites.xls", "xls", "PCoA样本坐标表", 0, ""],
            ["Nmds", "", "NMDS分析结果输出目录", 0, ""],
            ["Nmds/nmds_sites.xls", "xls", "NMDS样本坐标表", 0, ""],
            ["Nmds/nmds_stress.xls", "xls", "NMDS样本特征拟合度值", 0, ""],
            ["Hcluster", "", "Hcluster层级聚类分析结果目录", 0, ""],
            ["Hcluster/hcluster.tre", "tre", "Hcluster层级聚类树结果表", 0, ""],
            ["Hcluster/taxa.table.xls", "xls", "Hcluster各样本物种丰度结果表", 0, ""],
            ["Hcluster/taxa.precents.table.xls", "xls", "Hcluster各样本物种相对丰度结果表", 0, ""],
            ["Correlation", "", "Correlation相关性分析结果目录", 0, ""],
            ["Correlation/hcluster.tre", "tre", "Hcluster聚类树结果表", 0, ""],
        ]
        regexps = [
            [r"Correlation/correlation_*.xls", "xls", "Correlation相关性系数结果表", 0, ""],
            [r"Correlation/*_pvalue.xls", "xls", "Correlation相关性pvalues结果表", 0, ""]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(BetaDiversityWorkflow, self).end()
