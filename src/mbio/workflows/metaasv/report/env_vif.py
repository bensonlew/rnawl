# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os
import re
import json
import glob
from biocluster.workflow import Workflow
from mainapp.models.mongo.public.meta.meta import Meta
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir


class EnvVifWorkflow(Workflow):
    """
    metaasv VIF 方差膨胀因子分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnvVifWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "level", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "abund_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "asv_id", "type": "string"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "env_id", "type": "string"},
            {"name": "env_labs", "type": "string"},
            {"name": "viflim", "type": "int", "default": 10},
            {"name": "env_file", "type": "infile", 'format': "meta.otu.group_table"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.vif = self.add_tool('statistical.env_vif')
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")

    def run_sort_samples(self):
        """
        筛选丰度表
        :return:
        """
        self.logger.info("正常运行啦")
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
            self.sort_samples.set_options({
                "in_otu_table": abund_table,
                "group_table": self.option("group_file"),  ###
                "sample_del_warn": "T"
            })
            self.sort_samples.on("end", self.run_env_vif)
            self.sort_samples.run()

    def run_env_vif(self):
        """
        运行 env if
        :return:
        """
        abund_table = self.sort_samples.option("out_otu_table").prop['path']
        num_lines = open(abund_table, 'r').readlines()
        if len(num_lines) < 3:
            self.set_error('丰度表数据少于2行，请重新设置参数!')
        self.get_new_env_file()
        options = {
            "abundtable": self.sort_samples.option("out_otu_table"),
            "envtable": self.option('env_file'),
            "viflim": self.option("viflim")
        }
        self.vif.set_options(options)
        self.vif.on('end', self.set_db)
        self.vif.run()

    def set_db(self):
        """
        链接结果文件和导入MongoDB
        :return:
        """
        link_dir(self.vif.output_dir, self.output_dir)
        self.logger.info("正在写入mongo数据库")

        api_env_vif = self.api.api("metaasv.env_vif")
        env_result = os.listdir(self.output_dir)
        for file in env_result:
            if re.search("final_.*_vif.txt", file):
                final_path = os.path.join(self.output_dir, file)
                api_env_vif.add_env_vif_detail(final_path, "after", self.option("main_id"))
            elif re.search("raw_.*_vif.txt", file):
                raw_path = os.path.join(self.output_dir, file)
                api_env_vif.add_env_vif_detail(raw_path, "before", self.option("main_id"))
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        self.run_sort_samples()
        super(EnvVifWorkflow, self).run()

    def end(self):
        """
        结束和上传结果文件
        :return:
        """
        env_result = glob.glob(self.vif.output_dir + "/*")
        result_dir = self.add_upload_dir(self.output_dir)
        if re.search("cca_vif.txt", env_result[0]):
            result_dir.add_relpath_rules([
                [".", "", "VIF方差膨胀因子分析结果目录", 0, ""],
                ["./raw_cca_vif.txt", "txt", "筛选前VIF方差膨胀因子分析结果表", 0, ""],
                ["./final_cca_vif.txt", "txt", "筛选后VIF方差膨胀因子分析结果表", 0, ""],
                ["./DCA.txt", "txt", "判断计算VIF方差膨胀因子分析的方法文件", 0, ""]
            ])
        else:
            result_dir.add_relpath_rules([
                [".", "", "VIF方差膨胀因子分析结果目录", 0, ""],
                ["./raw_rda_vif.txt", "txt", "筛选前VIF方差膨胀因子分析结果表", 0, ""],
                ["./final_rda_vif.txt", "txt", "筛选后VIF方差膨胀因子分析结果表", 0, ""],
                ["./DCA.txt", "txt", "判断计算VIF方差膨胀因子分析的方法文件", 0, ""]
            ])
        super(EnvVifWorkflow, self).end()

    def get_new_env_file(self):
        """
        根据筛选的环境因子和环境因子表筛选环境因子
        :return:
        """
        file_path = os.path.join(self.work_dir, "selected_env_table.xls")
        envs = self.option('env_labs').split(',')
        fw = open(file_path, 'w')
        with open(self.option('env_file').path,'r') as f:
            heads = f.readline().strip().split('\t')
            env_index = [heads.index(i) for i in envs]
            get_index = [0]
            get_index.extend(env_index)
            fw.write('\t'.join([heads[i] for i in get_index])+'\n')
            for line in f:
                lines = line.strip().split('\t')
                fw.write('\t'.join([lines[i] for i in get_index])+'\n')
        self.option('env_file').prop['path'] = file_path

