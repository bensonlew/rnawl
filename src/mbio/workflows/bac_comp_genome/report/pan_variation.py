# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import os,shutil
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
from biocluster.core.exceptions import OptionError
import datetime
from biocluster.file import download, exists
from collections import defaultdict
import shutil
import gevent
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file


class PanVariationWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        细菌比较基因组泛基因组分析
        变异分析
        :return:
        """
        self._sheet = wsheet_object
        super(PanVariationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_dir", "type": "string"},  # 输入fasta的dir文件
            {"name": "cluster_path", "type": "infile", "format": "sequence.profile_table"},  # 物种注释总览表
            {"name": "main_id", "type": "string"},#输入pan_category表的id
            {"name": "from_main_id", "type": "string"},#输入pan主表表的id
            {"name": "update_info", "type": "string"},
        ]
        self.panvariation = self.add_tool("bac_comp_genome.pan_variation")
        self.api_genome = self.api.api('bac_comp_genome.pan')
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
        if not self.option("from_main_id"):
            raise OptionError("必须设置输入的pan_id")
        if not self.option("cluster_path").is_set:
            raise OptionError("请提供输入的聚类二维表！")
        else:
            self.logger.info(self.option('cluster_path').prop['path'])
        if not self.option("gene_dir"):
            raise OptionError("请提供输入序列的文件夹！")
        else:
            self.logger.info(self.option('gene_dir'))
        return True

    def run_variation(self):
        """
        对聚类结果进行变异分析
        :return:
        """
        self.logger.info("开始进行变异分析")
        self.panvariation.set_options({
            "cluster": self.option("cluster_path"),
            "infile_dir": self.fasta_dir,
            })
        self.panvariation.on("end", self.set_output)
        self.panvariation.run()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('设置结果目录')
        result_path = os.path.join(self.output_dir)
        link_dir(self.panvariation.output_dir,result_path)
        self.logger.info('设置结果目录成功')
        self.set_db()

    def set_db(self):
        """
        将数据导入mongo
        :return:
        """
        self.logger.info('正在写入mongo数据库')
        self.remote_dir = self._sheet.output
        api_pan = self.api.api('bac_comp_genome.pan_variation')
        pan_id = self.option("from_main_id")
        self.logger.info("正在导变异分析的结果表")
        main_id = self.option("main_id")
        variation_path = os.path.join(self.output_dir, 'CDS_variation_analysis.xls')
        var_number = 0
        if os.path.exists(variation_path):
            with open(variation_path, 'r') as var:
                for line in var:
                    var_number += 1
            if var_number > 1:
                api_pan.add_pan_variation_stat(main_id, variation_path,)
        self.end()

    def get_fasta_dir(self):
        """
        根据group和输入的文件夹获取此次分析的样本夹
        :return:
        """
        self.fasta_dir = os.path.join(self.work_dir, "fasta_dir")
        if os.path.exists(self.fasta_dir):
            shutil.rmtree(self.fasta_dir)
        os.mkdir(self.fasta_dir)
        specimen_list = self.api_genome.find_specimen(self.option("from_main_id"))
        for sample_name in specimen_list:
            for i in ["_CDS.fna", "_CDS.faa"]:
                if os.path.exists(self.option("gene_dir") + "/" + sample_name + "/" + sample_name + i):
                    self.logger.info("gene_dir: {}".format(self.option("gene_dir")))
                    if "_CDS.fnn" == i:
                        download(self.option("gene_dir") + "/" + sample_name + "/" + sample_name + i, self.fasta_dir+ "/" + sample_name + "_CDS.fna")
                    else:
                        download(self.option("gene_dir") + "/" + sample_name + "/" + sample_name + i, self.fasta_dir + "/" + sample_name + i)
                else:
                    download(self.option("gene_dir") + "/" + sample_name + "/" + sample_name + i, self.fasta_dir+ "/" + sample_name + i)

    def run(self):
        """
        开始运行了
        :return:
        """
        self.logger.info("开始运行")
        self.get_fasta_dir()
        self.run_variation()
        super(PanVariationWorkflow, self).run()

    def end(self):
        """
        结束了
        :return:
        """
        self.logger.info("开始结果文件上传")
        sdir = self.add_upload_dir(self.output_dir)
        repaths = [
            ["Pangenome", "", "同源基因突变分析输出目录",0,""],
            ["pangenomes/cluster/CDS_variation_analysis.xls", "", "同源基因突变分析结果表", 0],
            ["pangenomes/cluster/CDS_variation.xls", "", "同源基因突变分析统计表", 0],
            ]
        sdir.add_relpath_rules(repaths)
        super(PanVariationWorkflow, self).end()
