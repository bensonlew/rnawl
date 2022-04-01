# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

"""报告计算VIF方差膨胀因子"""

import os
import re
import json
import glob
from biocluster.workflow import Workflow
from mainapp.models.mongo.public.meta.meta import Meta
from biocluster.core.exceptions import OptionError
from mbio.packages.meta.common_function import envname_restore
from mbio.packages.meta.save_params import save_params


class EnvVifWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnvVifWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "level", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "params", "type": "string"},
           # {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},

            {"name": "abund_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "otu_id", "type": "string"},
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
        self.vif = self.add_tool('statistical.env_vif')
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")



    '''
    def run_get_abund_table(self):
        self.get_abund_table = self.add_tool("meta.create_abund_table")
        self.get_abund_table.set_options({
            'anno_table': self.option('anno_table'),
            'geneset_table': self.option('geneset_table'),
            'gene_list': self.option('gene_list'),
            'level_type': self.option('level_id'),
            'level_type_name': self.option('level_type_name'),
            'lowest_level': self.option('lowest_level'),
        })
        self.get_abund_table.on("end", self.run_sort_samples)
        self.get_abund_table.run()
    '''

    def run_sort_samples(self):
        self.logger.info("正常运行啦")
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            self.logger.info("缺少abund_file")
        #     self.logger.info("success")
        #    abund_table = self.get_abund_table.option("out_table").prop['path']
        self.logger.info(abund_table)
        # self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        #group_table_path = os.path.join(self.work_dir, "group_table.xls")  #guanqing.zou 20180419
        #group_table_path = Meta().group_detail_to_table(self.option("group_detail"), group_table_path) #guanqing.zou 20180419
        self.sort_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option("group_file"),  ###
            "sample_del_warn": "T"
        })
        self.sort_samples.on("end", self.run_env_vif)
        self.sort_samples.run()

    def run_env_vif(self):
        abund_table = self.sort_samples.option("out_otu_table").prop['path']
        num_lines = open(abund_table, 'r').readlines()
        if len(num_lines) < 3:
            self.set_error('丰度表数据少于2行，请重新设置参数!', code="12701201")
        self.get_new_env_file()  #guanqing.zou 20180425
        

        options = {
            "abundtable": self.sort_samples.option("out_otu_table"),
            "envtable": self.option('env_file'),
            "viflim": self.option("viflim"),
            
        }
        self.vif.set_options(options)
        self.vif.on('end', self.set_db)
        self.vif.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        api_env_vif = self.api.env_vif
        env_result = glob.glob(self.vif.output_dir + "/*")
        for file in env_result:
            if re.search("final_.*_vif.txt", file):
                final_path = file
                api_env_vif.add_env_vif_detail(final_path, "after_vif", self.option("main_id"))
            elif re.search("raw_.*_vif.txt", file):
                raw_path = file
                api_env_vif.add_env_vif_detail(raw_path, "before_vif", self.option("main_id"))
        self.end()

    def run(self):
        if self.option("abund_file").is_set:
            self.run_sort_samples()
        else:
            self.logger.info("缺少abund_file")
            #self.run_get_abund_table()
        self.output_dir = self.vif.output_dir
        super(EnvVifWorkflow, self).run()

    @envname_restore
    def end(self):
        env_result = glob.glob(self.vif.output_dir + "/*")
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        if re.search("cca_vif.txt", env_result[0]):
            result_dir.add_relpath_rules([
                [".", "", "VIF方差膨胀因子分析结果目录", 0, "110150"],
                ["./raw_cca_vif.txt", "txt", "筛选前VIF方差膨胀因子分析结果表", 0, "110152"],
                ["./final_cca_vif.txt", "txt", "筛选后VIF方差膨胀因子分析结果表", 0, "110151"],
                ["./DCA.txt", "txt", "判断计算VIF方差膨胀因子分析的方法文件", 0, "110153"]
            ])
        else:
            result_dir.add_relpath_rules([
                [".", "", "VIF方差膨胀因子分析结果目录", 0, "110150"],
                ["./raw_rda_vif.txt", "txt", "筛选前VIF方差膨胀因子分析结果表", 0, "110152"],
                ["./final_rda_vif.txt", "txt", "筛选后VIF方差膨胀因子分析结果表", 0, "110151"],
                ["./DCA.txt", "txt", "判断计算VIF方差膨胀因子分析的方法文件", 0, "110153"]
            ])
        super(EnvVifWorkflow, self).end()

    def get_new_env_file(self):
        ###guanqing.zou 20180425
        file_path = os.path.join(self.work_dir, "selected_env_table.xls")
        envs = self.option('env_labs').split(',')
        fw = open(file_path, 'w')
        with open(self.option('env_file').path,'r') as f:
        #self.logger.info('get_new_env_file'+self.option('env_file').path)
        #with open(os.path.join(self.work_dir, "env_file_input_env.xls"),'r') as f:
            heads = f.readline().strip().split('\t')
            env_index = [heads.index(i) for i in envs]
            get_index = [0]
            get_index.extend(env_index)
            fw.write('\t'.join([heads[i] for i in get_index])+'\n')
            for line in f:
                lines = line.strip().split('\t')
                fw.write('\t'.join([lines[i] for i in get_index])+'\n')
        self.option('env_file').prop['path'] = file_path
        #self.logger.info('get_new_env_file'+self.option('env_file').path)

