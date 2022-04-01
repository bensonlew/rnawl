# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import os
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir


class PanCategoryWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        细菌比较基因组泛基因组分析
        第二步 对聚类结果进行分组方案的合并
        :return:
        """
        self._sheet = wsheet_object
        super(PanCategoryWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "category", "type": "string"}, #输入分类方案类型{'core'：{"min":1, "max": 1}, 'dispensibale'：{"min":0.15, "max": 0.95}, 'soft_core'{"min":0.95, "max": 1}, 'unique'{"min":0, "max": 0.15}}min能取到，max取不到
            {"name": "percent", "type":"string", "default": "false"}, #false 表示按个数，true表示按百分比
            {"name": "category_name", "type":"string"}, #输入分组方案的名称
            {"name": "cluster_file", "type": "infile", "format": "sequence.profile_table"},  # 物种注释总览表
            # {"name": "group_id", "type": "string"},#输入group表的id
            {"name": "main_id", "type": "string"},#输入pan_category表的id
            {"name": "from_main_id", "type": "string"},#输入pan主表表的id
            # {"name": "class", "type": "string"}, #方案几
            {"name": "update_info", "type": "string"},
            {"name": "pan_group_id", "type": "string", 'default': ""}, # 分类方案的pan_group的id
            {"name": "pan_category_names", "type": "string", "default": ""}, # 分类方案的类型有哪些
        ]
        self.category = self.add_tool("bac_comp_genome.pan_category")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sample = []

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        self.logger.info("开始pan的参数检查")
        if not self.option("category"):
            raise OptionError("必须设置输入的category")
        if not self.option("cluster_file").is_set:
            raise OptionError("请提供输入序列的文件夹！")
        else:
            self.logger.info(self.option('cluster_file').prop['path'])
        return True

    def run_category(self):
        """
        根据分类方案计算整理
        :return:
        """
        self.logger.info("开始进行分组方案合并")
        opts = ({
            "cluster": self.option("cluster_file"),
            "category": self.option("category"),
            "category_name": self.option("category_name"),
            })
        if self.option("percent") in ["false", "False"]:
            opts["percent"] = False
        else:
            opts["percent"] = True
        self.category.set_options(opts)
        self.category.on("end", self.set_output)
        self.category.run()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('设置结果目录')
        result_path = os.path.join(self.output_dir)
        link_dir(self.category.output_dir,result_path)
        self.logger.info('设置结果目录成功')
        self.set_db()

    def set_db(self):
        """
        将数据导入mongo
        :return:
        """
        self.logger.info('正在写入mongo数据库')
        self.remote_dir = self._sheet.output
        api_category = self.api.api('bac_comp_genome.pan_category')
        pan_id = self.option("from_main_id")

        self.logger.info("正在导分组方案合并的主表")

        category_id = self.option("main_id")
        self.logger.info("正在导分组方案合并画图的详情表")
        distribution_path = os.path.join(self.output_dir, "pangenome_distribution.xls")
        api_category.add_category_graph(category_id, distribution_path)
        cluster_path = os.path.join(self.output_dir, "pangenome_clusters.xls")
        api_category.add_category_detail(category_id, cluster_path, pan_id=pan_id, type="cluster")

        gene_path = os.path.join(self.output_dir, "pangenome_genes.xls")
        api_category.add_category_detail(category_id, gene_path, pan_id=pan_id, type="gene")
        status_info = {"pan_group_id": self.option("pan_group_id"), "pan_category_names": self.option("pan_category_names")}
        api_category.update_sg_status(self.option('main_id'), status_info)
        self.end()


    def run(self):
        """
        开始运行了
        :return:
        """
        self.logger.info("开始运行")
        self.run_category()
        super(PanCategoryWorkflow, self).run()

    def end(self):
        """
        结束了
        :return:
        """
        self.logger.info("开始结果文件上传")
        sdir = self.add_upload_dir(self.output_dir)
        repaths = [
            ["Pangenome", "", "Pangenome分类结果输出目录",0,""],
            ["pangenome_clusters.xls", "xls", "Pangenome分类cluster统计表",0,""],
            ["pangenome_genes.xls", "xls", "Pangenome分类gene统计表",0,""],
            ["pangenome_distribution.xls", "xls", "Pangenome分类柱形图表",0,""],
            ]
        sdir.add_relpath_rules(repaths)
        super(PanCategoryWorkflow, self).end()
