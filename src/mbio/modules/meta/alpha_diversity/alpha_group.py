# -*- coding: utf-8 -*-


import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
from mainapp.models.mongo.public.meta.meta import Meta

class AlphaGroupModule(Module):

    def __init__(self, work_id):
        super(AlphaGroupModule, self).__init__(work_id)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的文件夹
            {"name": "estimate_indices", "type": "string", "default": "ace,chao,shannon,simpson,coverage"},
            {"name": "rarefy_indices", "type": "string", "default": "sobs,ace,chao,shannon,simpson,coverage"},  # 指数类型
            {"name": "rarefy_freq", "type": "int", "default": 100},
            {"name": "level", "type": "string", "default": "otu"},  # level水平
            {"name": "grouptable", "type": "string"},
        ]
        self.add_option(options)
        self.all_list = []
        self.all_list2 = []

    def run_all(self):
        for file in os.listdir(self.option("grouptable")):
            self.alpha = self.add_module('meta.alpha_diversity.alpha_diversity')
            options = {
                'otu_table': self.option('otu_table'),
                'level': self.option('level'),
                'estimate_indices': self.option('estimate_indices'),
                'rarefy_indices': self.option('rarefy_indices'),
                'rarefy_freq': self.option('rarefy_freq'),
                'group': self.option("grouptable") + "/" + file,
                'meta_group_name': file.rstrip(".group.txt")
            }
            self.alpha.set_options(options)
            self.all_list.append(self.alpha)
        if len(self.all_list) > 1:
            self.on_rely(self.all_list, self.run_pan_core)
        else:
            self.all_list[0].on('end', self.run_pan_core)
        for tool in self.all_list:
            tool.run()

    def run_pan_core(self):
        for file in os.listdir(self.option("grouptable")):
            self.pan_core = self.add_tool('meta.otu.pan_core_otu')
            options = {
                "in_otu_table": self.option('otu_table'),
                'group_table': self.option("grouptable") + "/" + file,
                'meta_group_name': file.rstrip(".group.txt")
            }
            self.pan_core.set_options(options)
            self.all_list2.append(self.pan_core)
        if len(self.all_list2) > 1:
            self.on_rely(self.all_list2, self.set_output)
        else:
            self.all_list2[0].on('end', self.set_output)
        for tool in self.all_list2:
            tool.run()

    def set_output(self):
        if self.all_list:
            for module in self.all_list:
                if not os.path.exists(self.output_dir + "/Alpha_diversity"):
                    os.mkdir(self.output_dir + "/Alpha_diversity")
                for dir in os.listdir(module.output_dir):
                    self.logger.info(module.output_dir + "/" + dir)
                    self.logger.info(self.output_dir + "/Alpha_diversity/")
                    self.logger.info(
                        'cp -r %s %s' % (module.output_dir + "/" + dir + "/", self.output_dir + "/Alpha_diversity/"))
                    if dir:
                        os.system('cp -r %s %s' % (
                        module.output_dir + "/" + dir + "/", self.output_dir + "/Alpha_diversity/"))
        if self.all_list2:
            for tool in self.all_list2:
                if not os.path.exists(self.output_dir + "/Pan_Core"):
                    os.mkdir(self.output_dir + "/Pan_Core")
                for dir in os.listdir(tool.output_dir):
                    if dir:
                        os.system('cp -r %s %s' % (tool.output_dir + "/" + dir + "/", self.output_dir + "/Pan_Core/"))
        self.end()

    def run(self):
        super(AlphaGroupModule, self).run()
        self.run_all()

    def end(self):
        super(AlphaGroupModule, self).end()

