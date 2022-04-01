# -*- coding: utf-8 -*-


import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
from mainapp.models.mongo.public.meta.meta import Meta

class DiffGroupModule(Module):

    def __init__(self, work_id):
        super(DiffGroupModule, self).__init__(work_id)
        options = [
            {"name": "otutable", "type": "string"},  # 输入的文件夹
            {"name": "grouptable", "type": "string"}, # group 文件的目录

            {"name": "multi_correction", "type": "string", "default": "fdr"},
            {"name": "multi_test", "type": "string", "default": "kru_H"},
            {"name": "multi_methor", "type": "string", "default": "tukeykramer"},
            {"name": "multi_coverage", "type": "float", "default": 0.95},

            {"name": "two_ci", "type": "float", "default": 0.05},
            {"name": "two_correction", "type": "string", "default": "fdr"},
            {"name": "two_test", "type": "string", "default": "mann"},
            {"name": "two_type", "type": "string", "default": "two.side"},
            {"name": "two_coverage", "type": "float", "default": 0.95},

            {"name": "lefse_type", "type": "string", "default": "meta_taxon"},
            {"name": "lda_filter", "type": "float", "default": 2.0},
            {"name": "strict", "type": "int", "default": 1},
            {"name": "start_level", "type": "int", "default": 3},
            {"name": "end_level", "type": "int", "default": 7},
            {"name": "normalization", "type": "int", "default": 0},  # 是否标准化 1，0
        ]
        self.add_option(options)
        self.multi_phy = []
        self.multi_gen = []
        self.two_phy = []
        self.phy = []
        self.genus = []
        self.two_gen = []
        self.lefse_group = []

    def run_phy(self):
        for file in os.listdir(self.option("grouptable")):
            group_name = []
            with open(self.option("grouptable")+"/"+file) as f:
                data = f.readlines()
                group_gname = data[0].strip().split("\t")[1]
                for i in data[1:]:
                    group_name.append(i.strip().split("\t")[1])

            if len(set(group_name)) == 2:
                self.two_group = self.add_tool("statistical.metastat")
                options = {
                    "mann_input": self.option("otutable")+"/otu_taxon_Phylum.xls",
                    "mann_group": self.option("grouptable") + "/" + file,
                    "mann_ci": self.option("two_ci"),
                    "mann_correction": self.option("two_correction"),
                    "mann_type": self.option("two_type"),
                    "test": self.option("two_test"),
                    "mann_gname": group_gname,
                    "mann_coverage": self.option("two_coverage"),
                    "meta_group_name": file.rstrip(".group.txt") + "/" + "DiffStatTwoGroup_Phylum"
                }
                self.two_group.set_options(options)
                self.phy.append(self.two_group)

            elif len(set(group_name)) >2:
                self.multiple = self.add_tool("statistical.metastat")
                options = {
                    "kru_H_input": self.option("otutable")+"/otu_taxon_Phylum.xls",
                    "kru_H_group": self.option("grouptable") + "/" + file,
                    "kru_H_correction": self.option("multi_correction"),
                    "test": self.option("multi_test"),
                    "kru_H_gname": group_gname,
                    "kru_H_methor": self.option("multi_methor"),
                    "kru_H_coverage": self.option("multi_coverage"),
                    "meta_group_name": file.rstrip(".group.txt") + "/" + "DiffStatMultiple_Phylum"
                }
                self.multiple.set_options(options)
                self.phy.append(self.multiple)
            else:
                pass
        if len(self.phy) > 1:
            self.on_rely(self.phy, self.run_genu)
        elif len(self.phy) == 0:
            self.run_lefse()
        else:
            self.phy[0].on('end', self.run_genu)
        for tool in self.phy:
            tool.run()

    def run_genu(self):
        for file in os.listdir(self.option("grouptable")):
            group_name = []
            with open(self.option("grouptable") + "/" + file) as f:
                data = f.readlines()
                group_gname = data[0].strip().split("\t")[1]
                for i in data[1:]:
                    group_name.append(i.strip().split("\t")[1])

            if len(set(group_name)) == 2:
                self.two_group = self.add_tool("statistical.metastat")
                options = {
                    "mann_input": self.option("otutable")+"/otu_taxon_Genus.xls",
                    "mann_group": self.option("grouptable") + "/" + file,
                    "mann_ci": self.option("two_ci"),
                    "mann_correction": self.option("two_correction"),
                    "mann_type": self.option("two_type"),
                    "test": self.option("two_test"),
                    "mann_gname": group_gname,
                    "mann_coverage": self.option("two_coverage"),
                    "meta_group_name": file.rstrip(".group.txt") + "/" + "DiffStatTwoGroup_Genus"
                }
                self.two_group.set_options(options)
                self.genus.append(self.two_group)

            elif len(set(group_name)) > 2:
                self.multiple = self.add_tool("statistical.metastat")
                options = {
                    "kru_H_input": self.option("otutable")+"/otu_taxon_Genus.xls",
                    "kru_H_group": self.option("grouptable") + "/" + file,
                    "kru_H_correction": self.option("multi_correction"),
                    "test": self.option("multi_test"),
                    "kru_H_gname": group_gname,
                    "kru_H_methor": self.option("multi_methor"),
                    "kru_H_coverage": self.option("multi_coverage"),
                    "meta_group_name": file.rstrip(".group.txt") + "/" + "DiffStatMultiple_Genus"
                }
                self.multiple.set_options(options)
                self.genus.append(self.multiple)
            else:
                pass
        if len(self.genus) > 1:
            self.on_rely(self.genus, self.run_lefse)
        elif len(self.genus) == 0:
            self.run_lefse()
        else:
            self.genus[0].on('end', self.run_lefse)
        for tool in self.genus:
            tool.run()

    def run_lefse(self):
        for file in os.listdir(self.option("grouptable")):
            group_name = []
            with open(self.option("grouptable") + "/" + file) as f:
                data = f.readlines()
                group_gname = data[0].strip().split("\t")[1]
                for i in data[1:]:
                    group_name.append(i.strip().split("\t")[1])
            self.lefse = self.add_tool("statistical.lefse")
            options = {
                "lefse_type": self.option("lefse_type"),
                "lefse_input": self.option("otutable").strip("tax_summary_a")+"otu_summary.xls",
                "lefse_group": self.option("grouptable") + "/" + file,
                "lda_filter": self.option("lda_filter"),
                "strict": self.option("strict"),
                "lefse_gname": group_gname,
                "start_level": self.option("start_level"),
                "end_level": self.option("end_level"),
                "percent_abund": "true",
                "meta_group_name": file.rstrip(".group.txt")
            }
            if not self.option('normalization'):
                options['percent_abund'] = 'false'
            self.lefse.set_options(options)
            self.lefse_group.append(self.lefse)
        if len(self.lefse_group) > 1:
            self.on_rely(self.lefse_group, self.set_output)
        else:
            self.lefse_group[0].on('end', self.set_output)
        for tool in self.lefse_group:
            tool.run()

    def set_output(self):
        for module in self.phy:
            os.system('cp -r %s %s' % (module.output_dir, self.output_dir))
        for module in self.genus:
            os.system('cp -r %s %s' % (module.output_dir, self.output_dir))
        for module in self.lefse_group:
            os.system('cp -r %s %s' % (module.output_dir, self.output_dir))
        self.end()

    def run(self):
        super(DiffGroupModule, self).run()
        self.run_phy()

    def end(self):
        super(DiffGroupModule, self).end()
