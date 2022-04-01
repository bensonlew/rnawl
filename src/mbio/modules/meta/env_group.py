# -*- coding: utf-8 -*-


import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
from mainapp.models.mongo.public.meta.meta import Meta

class EnvGroupModule(Module):

    def __init__(self, work_id):
        super(EnvGroupModule, self).__init__(work_id)
        options = [
            {"name": "otutable", "type": "string"},  # 输入的文件夹
            {"name": "grouptable", "type": "string"}, # group 文件的目录
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
        ]
        self.add_option(options)
        self.rdacca_list = []
        self.heatmap_list = []
        self.heatmap_list_genus = []
        self.mantel_list = []


    def run_rdacca(self):
        for file in os.listdir(self.option("grouptable")):
            sample_num = len(open(self.option("grouptable") + "/" + file).readlines())
            if sample_num > self.env_len:
                self.rdacca = self.add_tool('meta.beta_diversity.rda_cca')
                options = {
                    "otutable": self.option("otutable") + "/otu_taxon_otu.full.xls",
                    "envtable": self.option("envtable"),
                    "group_table": self.option("grouptable") + "/" + file,
                    "meta_group_name": file.rstrip(".group.txt")
                }
                self.rdacca.set_options(options)
                self.rdacca_list.append(self.rdacca)
        if len(self.rdacca_list) > 1:
            self.on_rely(self.rdacca_list, self.run_heatmap_phy)
        elif len(self.rdacca_list) == 0:
            self.run_heatmap_phy()
        else:
            self.rdacca_list[0].on('end', self.run_heatmap_phy)
        for tool in self.rdacca_list:
            tool.run()

    def run_heatmap_phy(self):
        for file in os.listdir(self.option("grouptable")):
            if self.env_len >= 2:
                self.correlation = self.add_tool('statistical.pearsons_correlation')
                options = {
                    "otutable": self.option("otutable") + "/otu_taxon_Phylum.xls",
                    "envtable": self.option("envtable"),
                    "env_cluster": "average",
                    "species_cluster": "average",
                    "top_species": 50,
                    "method": "spearmanr",
                    "level": "phylum",
                    "meta_group_name": file.rstrip(".group.txt")+"/SpearmanCorrelation_Phylum"
                }
                self.correlation.set_options(options)
                self.heatmap_list.append(self.correlation)
        if len(self.heatmap_list) > 1:
            self.on_rely(self.heatmap_list, self.run_heatmap_genus)
        elif len(self.heatmap_list) == 0:
            self.run_heatmap_genus()
        else:
            self.heatmap_list[0].on('end', self.run_heatmap_genus)
        for tool in self.heatmap_list:
            tool.run()

    def run_heatmap_genus(self):
        for file in os.listdir(self.option("grouptable")):
            if self.env_len >= 2:
                self.correlation = self.add_tool('statistical.pearsons_correlation')
                options = {
                    "otutable": self.option("otutable")+ "/otu_taxon_Genus.xls",
                    "envtable": self.option("envtable"),
                    "env_cluster": "average",
                    "species_cluster": "average",
                    "top_species": 50,
                    "level": "genus",
                    "meta_group_name": file.rstrip(".group.txt") + "/SpearmanCorrelation_Genus"
                }
                self.correlation.set_options(options)
                self.heatmap_list_genus.append(self.correlation)
        if len(self.heatmap_list_genus) > 1:
            self.on_rely(self.heatmap_list_genus, self.run_mantel)
        elif len(self.heatmap_list_genus) == 0:
            self.run_mantel()
        else:
            self.heatmap_list_genus[0].on('end', self.run_mantel)
        for tool in self.heatmap_list_genus:
            tool.run()

    def run_mantel(self):
        for file in os.listdir(self.option("grouptable")):
            self.mantel = self.add_module('statistical.mantel_test')
            options = {
                "otutable": self.option("otutable") + "/otu_taxon_otu.full.xls",
                "factor": self.option("envtable").prop["path"],
                "otumatrixtype": "bray_curtis",
                "factormatrixtype": "bray_curtis",
                "meta_group_name": file.rstrip(".group.txt") + "/MantelTest"
            }
            self.mantel.set_options(options)
            self.mantel_list.append(self.mantel)
        if len(self.mantel_list) > 1:
            self.on_rely(self.mantel_list, self.set_output)
        else:
            self.mantel_list[0].on('end', self.set_output)
        for tool in self.mantel_list:
            tool.run()

    def set_output(self):
        if self.rdacca_list:
            for module in self.rdacca_list:
                for dir in os.listdir(module.output_dir):
                    if dir:
                        os.system('cp -r %s %s' % (module.output_dir, self.output_dir))
        if self.heatmap_list:
            for module in self.heatmap_list:
                for dir in os.listdir(module.output_dir):
                    if dir:
                        os.system('cp -r %s %s' % (module.output_dir, self.output_dir))
        if self.heatmap_list_genus:
            for module in self.heatmap_list_genus:
                for dir in os.listdir(module.output_dir):
                    if dir:
                        os.system('cp -r %s %s' % (module.output_dir, self.output_dir))
        if self.mantel_list:
            for module in self.mantel_list:
                for dir in os.listdir(module.output_dir):
                    if dir:
                        os.system('cp -r %s %s' % (module.output_dir, self.output_dir))
        self.end()

    def run(self):
        super(EnvGroupModule, self).run()
        self.env_len = 0
        with open(self.option("envtable").path) as v:
            self.env_len = len(v.readline().strip().split("\t")[1:])
        self.run_rdacca()


    def end(self):
        super(EnvGroupModule, self).end()
