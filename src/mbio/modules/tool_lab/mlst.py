# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2021.04.08

import os,shutil,re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file
import pandas as pd

class MlstModule(Module):
    """
    多个样品的MLST进化树
    """

    def __init__(self, work_id):
        super(MlstModule, self).__init__(work_id)
        option = [
            {"name": "fasta_dir", "type": "infile", "format": "tool_lab.fasta_dir"},  # 基因组文件夹
            {"name": "species", "type": "string", },  # 物种名称
        ]
        self.add_option(option)
        self.mlst_tree =self.add_tool("tool_lab.mlst_tree")
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("fasta_dir").is_set:
            raise OptionError("必须提fasta_dir！")

        if not self.option("species"):
            raise OptionError("必须提供species！")
        self.samples,self.fastas = self.get_sample(self.option("fasta_dir").prop['path']+"/list.txt")
        if len(self.samples) <5:
            raise OptionError("MLST进化树分析需要5个样本以上！")
        return True

    def run_mlst(self):
        for sample in self.samples:
            mlst = self.add_tool('tool_lab.mlst')
            opts = {
                "fasta": self.option("fasta_dir").prop['path'] + "/" + self.fastas[sample],
                "species": self.option("species"),
                "sample_name": sample,
            }
            mlst.set_options(opts)
            self.all_modules.append(mlst)
        if len(self.all_modules) > 1:
            self.on_rely(self.all_modules, self.run_tree)
        for module in self.all_modules:
            module.run()

    def run_tree(self):
        if not os.path.exists(self.work_dir + "/tmp_mlst"):
            os.mkdir(self.work_dir + "/tmp_mlst")
        if len(self.all_modules) > 1:
            for module in self.all_modules:
                for file in os.listdir(module.output_dir):
                    if re.search("mlst.ST.xls",file):
                        link_file(module.output_dir+"/"+file, self.work_dir + "/tmp_mlst/"+file)
        list2 = []
        for file in os.listdir(self.work_dir + "/tmp_mlst"):
            gene = pd.read_table(self.work_dir + "/tmp_mlst/" + file, sep='\t', header=0)
            list2.append(gene)
        all_data = pd.concat(list2)
        funs2 = all_data.columns.tolist()
        funs2.remove('sample')
        funs2.remove('ST')
        samples1 = ['sample', "ST"] + funs2
        all_data.to_csv(self.work_dir + "/all.ST.xls", sep='\t', columns=samples1, header=True, index=False)
        del all_data["ST"]
        funs = all_data.columns.tolist()
        funs.remove('sample')
        def cent_fun(values):
            num =0
            for v in values:
                if v == "-":
                    num +=1
            if num+1 == len(values):
                return 0
            else:
                return 1
        all_data = all_data[all_data.apply(lambda x: cent_fun(x) == 1, axis=1)]
        if all_data.shape[0] >=5:
            samples =['sample'] + funs
            all_data.to_csv(self.work_dir + "/all.martix.xls", sep='\t', columns=samples, header=True, index=False)
            self.mlst_tree.set_options({
                "martix": self.work_dir + "/all.martix.xls",
            })
            self.mlst_tree.on("end", self.set_output)
            self.mlst_tree.run()
        else:
            self.set_output()


    def set_output(self):
        link_file(self.work_dir + "/all.ST.xls",self.output_dir+"/all.ST.xls")
        if os.path.exists(self.mlst_tree.output_dir+"/mlst_tree.nwk"):
            link_file(self.mlst_tree.output_dir+"/mlst_tree.nwk", self.output_dir+"/mlst_tree.nwk")
        self.end()

    def get_sample(self,file):
        list = []
        dict = {}
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list.append(lin[0])
                dict[lin[0]] = lin[1]
        return list,dict

    def run(self):
        super(MlstModule, self).run()
        self.run_mlst()

    def end(self):
        super(MlstModule, self).end()