# -*- coding: utf-8 -*-


import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
from mainapp.models.mongo.public.meta.meta import Meta

class CompositionGroupModule(Module):

    def __init__(self, work_id):
        super(CompositionGroupModule, self).__init__(work_id)
        options = [
            {"name": "otutable", "type": "string"},  # 输入的文件夹
            {"name": "grouptable", "type": "string"}, # group 文件的目录
        ]
        self.add_option(options)
        self.comp_phy = []
        self.comp_genus = []
        self.comp_venn = []

    def run_phylum(self):
        for file in os.listdir(self.option("grouptable")):
            self.comp_analysis = self.add_module('meta.composition.composition_analysis_meta')
            options = {
                'analysis': "bar,circos",
                'abundtable': self.option("otutable")+"/otu_taxon_Phylum.full.xls",
                'group': self.option("grouptable") + "/" + file,
                'sample_method': "average",
                'method': "average",
                'meta_group_name': file.rstrip(".group.txt")
            }
            self.comp_analysis.set_options(options)
            self.comp_phy.append(self.comp_analysis)
        if len(self.comp_phy) > 1:
            self.on_rely(self.comp_phy, self.run_genus)
        elif len(self.comp_phy) == 0:
            self.run_genus()
        else:
            self.comp_phy[0].on('end', self.run_genus)
        for tool in self.comp_phy:
            tool.run()

    def run_genus(self):
        for file in os.listdir(self.option("grouptable")):
            self.comp_analysis = self.add_module('meta.composition.composition_analysis_meta')
            options = {
                'analysis': "bar,heatmap",
                'abundtable': self.option("otutable")+"/otu_taxon_Genus.full.xls",
                'group': self.option("grouptable") + "/" + file,
                'meta_group_name': file.rstrip(".group.txt"),
                'sample_method': "average",
                'method': "average",
            }
            self.comp_analysis.set_options(options)
            self.comp_genus.append(self.comp_analysis)
        if len(self.comp_genus) > 1:
            self.on_rely(self.comp_genus, self.run_venn)
        elif len(self.comp_genus) == 0:
            self.run_venn()
        else:
            self.comp_genus[0].on('end', self.run_venn)
        for tool in self.comp_genus:
            tool.run()

    def run_venn(self):
        for file in os.listdir(self.option("grouptable")):
            all_sample = []
            group_dict = {}
            with open(self.option("grouptable") + "/" + file) as f:
                for i in f.readlines()[1:]:
                    all_sample.append(i.strip().split("\t")[0])
                    if i.strip().split("\t")[1] in group_dict.keys():
                        group_dict[i.strip().split("\t")[1]].append(i.strip().split("\t")[0])
                    else:
                        group_dict[i.strip().split("\t")[1]] = [i.strip().split("\t")[0]]
            with open(self.option("otutable") + "/otu_taxon_otu.full.xls") as f1,open(self.work_dir+"/abundtable_tmp.txt","w") as t:
                data = f1.readlines()
                t.write(data[0])
                for i in data[1:]:
                    ii = i.strip().replace("; ",";")
                    iii = ii.replace(";","; ")
                    t.write(iii+"\n")
            if len(group_dict) >= 2:
                self.comp_analysis = self.add_module('meta.composition.composition_analysis_meta')
                options = {
                    'analysis': "venn",
                    'abundtable': self.work_dir+"/abundtable_tmp.txt",
                    'group': self.option("grouptable") + "/" + file,
                    'meta_group_name': file.rstrip(".group.txt")
                }
                self.comp_analysis.set_options(options)
                self.comp_venn.append(self.comp_analysis)
        if len(self.comp_venn) > 1:
            self.on_rely(self.comp_venn, self.set_output)
        elif len(self.comp_venn) == 0:
            self.set_output()
        else:
            self.comp_venn[0].on('end', self.set_output)
        for tool in self.comp_venn:
            tool.run()

    def set_output(self):
        for module in self.comp_phy:
            for dir in os.listdir(module.output_dir):
                if dir:
                    os.system('cp -r %s %s' % (module.output_dir+"/"+dir, self.output_dir))
        for module in self.comp_genus:
            for dir in os.listdir(module.output_dir):
                if dir:
                    os.system('cp -r %s %s' % (module.output_dir+"/"+dir, self.output_dir))
        for module in self.comp_venn:
            for dir in os.listdir(module.output_dir):
                if dir:
                    os.system('cp -r %s %s' % (module.output_dir+"/"+dir, self.output_dir))
        self.end()

    def run(self):
        super(CompositionGroupModule, self).run()
        self.run_phylum()

    def end(self):
        super(CompositionGroupModule, self).end()

