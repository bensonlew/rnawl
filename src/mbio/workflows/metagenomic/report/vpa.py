# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError
from comm_table import CommTableWorkflow
from mbio.packages.meta.common_function import envname_restore
import json
import collections

class VpaWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VpaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_type", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name":"params","type":"string"},
            {"name":"update_info","type":"string"},
            {"name":"main_id","type":"string"},
            {"name":"group","type":"infile","format": "meta.otu.group_table"},
            {"name": "env_file", "type": "infile", "format": "meta.otu.otu_table"} , #guanqing.zou 20180926
            {"name": "env_group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_detail", "type": "string"},
            {"name": "env_id", "type": "string"}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        #self.get_func_abund_table = self.add_tool('meta.create_abund_table')
        self.sort_func_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.vpa = self.add_tool('statistical.vpa')
        self.common_api = self.api.api("metagenomic.common_api")

    def run_func_sort_samples(self):
        abund_table = self.abundance.option("out_table").prop['path']
        self.sort_func_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group'),
        })
        self.sort_func_samples.on("end", self.run_vpa)
        self.sort_func_samples.run()

    def rewrite_env_group(self):
        self.new_env_group = os.path.join(self.work_dir, "env_group_new.txt")
        print self.option("params")
        env_group_id = json.loads(self.option("params"))['env_group_id']
        env_group = self.common_api.db['env_group'].find_one({"_id": ObjectId(env_group_id)})
        group_names = env_group["category_names"]
        env_names = env_group["specimen_names"]
        env_group_dict = collections.OrderedDict()
        for i in range(len(group_names)):
            env_group_dict[group_names[i]] = env_names[i]
        env_id_list = []
        with open(self.option("env_group").prop['path'], "r") as f:
            data = f.readlines()
            for lines in data[1:]:
                line = lines.strip().split("\t")
                env_id = line[0]
                env_id_list.append(env_id)
        with open(self.new_env_group, "w") as w:
            w.write("#env\tgroup\n")
            for i in env_group_dict.keys():
                for j in env_group_dict[i]:
                    if j in env_id_list:
                        w.write("{}\t{}\n".format(j, i))

    def run_vpa(self):
        func_abund_table = self.sort_func_samples.option("out_otu_table")
        self.vpa.set_options({
            "species_table" :  func_abund_table,
            "env_table" : self.option("env_file"),
            "group_table": self.new_env_group
            # "group_table": self.option('env_group')
        })
        self.vpa.on("end", self.set_db)
        self.vpa.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        #self.run_get_func_abund_table()
        self.rewrite_env_group()
        self.run_abundance(self.run_func_sort_samples)
        self.abundance.run()
        super(VpaWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        self.api_vpa = self.api.api("metagenomic.vpa")

        data_path = self.vpa.output_dir + '/env.R2adj.xls'
        graph_path = self.vpa.work_dir + '/env.plot.xls'
        self.api_vpa.add_vpa_detail(data_path, self.option("main_id"))
        png_path = self._sheet.output
        self.api_vpa.add_vpa_graph(graph_path, self.option("main_id"), png_path)
        self.output_dir = self.vpa.output_dir
        self.end()

    @envname_restore
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "VPA结果文件目录",0,"120341"],
            ["./env.R2adj.xls", "xls", "R2adj数据",0,"120342"],
            #["./env.plot.xls", "xls", "画图数据", 0],
            ["./vpa.png", "png", "vpa图",0],
            ["./vpa.pdf", "pdf", "vpa图",0]
        ])
        super(VpaWorkflow, self).end()
