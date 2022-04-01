# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last_modify:20180906

import os
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class TreeGenericModule(Module):
    """
    进化树--maximum_likehood(ML), neighbor_joining(NJ), Bayes
    """
    def __init__(self, work_id):
        super(TreeGenericModule, self).__init__(work_id)
        options = [
            {"name": "recode_vcf_path", "type": "infile", "format": "dna_gmap.vcf", "required": True},
            {"name": "tree_type", "type": "string", "default": "nj"},
            {"name": "bs_trees", "type": "int", "default": 1000}
        ]
        self.add_option(options)
        self.vcf2tree = self.add_tool("dna_evolution.vcf2tree")
        self.modeltest = self.add_tool("dna_evolution.modeltest")
        self.ml_tree = None
        self.nj_tree = None
        self.bayes_tree = None

    def check_options(self):
        if self.option("tree_type").lower() not in ["nj", 'ml', 'bayes']:
            raise OptionError("树的计算方法必须是nj or ml or bayes", code="11111")
        if type(self.option("bs_trees")) != int:
            raise OptionError("bs_trees参数必须是整型", code="11111")
        else:
            if self.option("bs_trees") <= 1:
                raise OptionError("bs_trees参数必须是大于1的整数", code="11111")

    def vcf2tree_run(self):
        self.vcf2tree.set_options({
            "recode_vcf_path": self.option("recode_vcf_path").prop['path']
        })
        self.vcf2tree.on("end", self.set_output, "vcf2tree")
        self.vcf2tree.run()

    def modeltest_run(self):
        self.modeltest.set_options({
            "phylip": self.vcf2tree.option("phylip_path")
        })
        self.modeltest.on("end", self.set_output, "modeltest")
        self.modeltest.run()

    def nj_tree_run(self):
        self.nj_tree.set_options({
            "pop_fasta": self.vcf2tree.option("pop_fasta")
        })
        self.nj_tree.on("end", self.set_output, "nj_tree")
        self.nj_tree.run()

    def ml_tree_run(self):
        self.ml_tree.set_options({
            "model_test_out": self.modeltest.option("model_test"),
            "bs_trees": self.option("bs_trees")
        })
        self.ml_tree.on("end", self.set_output, "ml_tree")
        self.ml_tree.run()

    def bayes_tree_run(self):
        self.bayes_tree.set_options({
            "pop_fasta": self.vcf2tree.option("pop_fasta"),
            "model_test_out": self.modeltest.option("model_test")
        })
        self.bayes_tree.on("end", self.set_output, "bayes_tree")
        self.bayes_tree.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'vcf2tree':
            self.linkdir(obj.output_dir, 'vcf2tree')
        elif event['data'] == 'modeltest':
            self.linkdir(obj.output_dir, 'modeltest')
        elif event['data'] == 'nj_tree':
            self.linkdir(obj.output_dir, 'tree')
        elif event['data'] == 'ml_tree':
            self.linkdir(obj.output_dir, 'tree')
        elif event['data'] == 'bayes_tree':
            self.linkdir(obj.output_dir, 'tree')
        else:
            pass

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        super(TreeGenericModule, self).run()
        if self.option("tree_type").lower() == "ml":
            self.ml_tree = self.add_tool("dna_evolution.ml_tree")
            self.vcf2tree.on("end", self.modeltest_run)
            self.modeltest.on("end", self.ml_tree_run)
            self.ml_tree.on("end", self.end)
        elif self.option("tree_type").lower() == "bayes":
            self.bayes_tree = self.add_tool("dna_evolution.bayes_tree")
            self.vcf2tree.on("end", self.modeltest_run)
            self.modeltest.on("end", self.bayes_tree_run)
            self.bayes_tree.on("end", self.end)
        else:
            self.nj_tree = self.add_tool("dna_evolution.nj_tree")
            self.vcf2tree.on("end", self.nj_tree_run)
            self.nj_tree.on("end", self.end)
        self.vcf2tree_run()

    def end(self):
        super(TreeGenericModule, self).end()
