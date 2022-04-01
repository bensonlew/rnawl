# -*- coding: utf-8 -*-


import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
from mainapp.models.mongo.public.meta.meta import Meta
import pandas as pd

class BetaGroupModule(Module):

    def __init__(self, work_id):
        super(BetaGroupModule, self).__init__(work_id)
        options = [
            {"name": "analysis", "type": "string",
             "default": "distance,anosim,pca,pcoa,nmds,rda_cca,dbrda,hcluster,plsda"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            {"name": "dbrda_method", "type": "string", "default": ""},
            # 当设定此值时，dbrda的计算方式将会改变，使用R中自带的距离算法，而不是先计算好距离矩阵，此处的计算方式与一般的距离计算的的值不一致
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "permutations", "type": "int", "default": 999},
            {"name": "linkage", "type": "string", "default": "average"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "envlabs", "type": "string", "default": ""},
            {"name": "pca_envlabs", "type": "string", "default": ""},
            {"name": "dbrda_envlabs", "type": "string", "default": ""},
            {"name": "rda_envlabs", "type": "string", "default": ""},
            {"name": "grouptable", "type": "string"},
            {"name": "grouplab", "type": "string", "default": ""},
            {"name": "anosim_grouplab", "type": "string", "default": ""},
            {"name": "plsda_grouplab", "type": "string", "default": ""},
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"},
            {"name": "dis_newicktree", "type": "outfile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "scale", "type": "string", "default": "F"},  # pca是否进行标准化 ，add by zouxuan
            {"name": "ellipse", "type": "string", "default": "F"},
            {"name": "diff_test_method", "type": "string", "default": ""},  # by houshuang 20190924 pca/pcoa/nmds组间差异检验
            {"name": "change_times", "type": "string", "default": ""},  # by houshuang 20190924 pca/pcoa/nmds组间差异检验
            {"name": "others_value", "type": "float", "default": "0.05"}
        ]
        self.add_option(options)
        self.all_list = []

    def run_all(self):
        num = 1
        for file in os.listdir(self.option("grouptable")):
            self.beta = self.add_module('meta.beta_diversity.beta_diversity')
            analysis = self.option('analysis')
            all_sample = []
            group_dict = {}
            with open(self.option("grouptable") + "/" + file) as f:
                for i in f.readlines()[1:]:
                    all_sample.append(i.strip().split("\t")[0])
                    if i.strip().split("\t")[1] in group_dict.keys():
                        group_dict[i.strip().split("\t")[1]].append(i.strip().split("\t")[0])
                    else:
                        group_dict[i.strip().split("\t")[1]] = [i.strip().split("\t")[0]]
            if len(group_dict.keys()) >=2:
                for x in group_dict.keys():
                    if len(group_dict[x]) ==1:
                        analysis = analysis.replace(",anosim","")
            else:
                analysis = analysis.replace(",anosim", "")
            if len(all_sample) < 3:
                analysis = "distance,hcluster"
            self.logger.info('analysis:{}'.format(analysis))
            df = pd.read_table(self.option('otutable').prop["path"],sep='\t',header=0)
            raw_sample = df.columns.values
            for xx in range(len(raw_sample)-1):
                if raw_sample[xx+1] not in all_sample:
                    df.drop(raw_sample[xx+1], axis=1, inplace=True)
            df.to_csv(self.work_dir+ "/otutable"+str(num),sep="\t",index=0)
            options = {
                'analysis': analysis,
                'dis_method': self.option('dis_method'),
                'otutable': self.work_dir+ "/otutable"+str(num),
                'permutations': self.option('permutations'),
                'group': self.option("grouptable") + "/" + file,
                'meta_group_name': file.rstrip(".group.txt"),
                'diff_test_method': self.option("diff_test_method"),
                'others_value': self.option("others_value")
            }
            self.beta.set_options(options)
            self.all_list.append(self.beta)
            num += 1
        if len(self.all_list) > 1:
            self.on_rely(self.all_list, self.set_output)
        else:
            self.all_list[0].on('end', self.set_output)
        for tool in self.all_list:
            tool.run()

    def set_output(self):
        if self.all_list:
            if not os.path.exists(self.output_dir + "/Beta_diversity/"):
                os.mkdir(self.output_dir + "/Beta_diversity/")
            for module in self.all_list:
                for dir in os.listdir(module.output_dir):
                    if dir:
                        os.system('cp -r %s %s' % (module.output_dir + "/" + dir, self.output_dir + "/Beta_diversity/"))
        self.end()

    def run(self):
        super(BetaGroupModule, self).run()
        self.run_all()

    def end(self):
        super(BetaGroupModule, self).end()

