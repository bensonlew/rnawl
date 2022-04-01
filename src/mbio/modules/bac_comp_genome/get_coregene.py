# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.10.09

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file,add_merge,add_fasta

class GetCoregeneModule(Module):
    """
    多个样品的core gene的获取
    """

    def __init__(self, work_id):
        super(GetCoregeneModule, self).__init__(work_id)
        option = [
            {"name": "gene_path", "type": "string"},  # string类型的基因预测dir路径
            {"name": "sample_list", "type": "infile", "format": "sequence.profile_table"},#文件中有样品信息
        ]
        self.add_option(option)
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("sample_list").is_set:
            raise OptionError("必须提供样品的信息文件！")
        if not self.option("gene_path"):
            raise OptionError("必须提供gene_path的路径！")
        return True

    def run_all_coregene(self):
        self.samples = self.get_sample(self.option("sample_list").prop['path'])
        for sample in self.samples:
            core_gene = self.add_tool('bac_comp_genome.get_core_gene')
            opts = {
                "seq_faa": self.option("gene_path") + sample + "/" + sample + "_CDS.faa",
                "seq_gff": self.option("gene_path") + sample + "/" + sample + "_CDS.gff",
                "sample_name": sample,
                "method": "genome",
            }
            core_gene.set_options(opts)
            self.all_modules.append(core_gene)
        if len(self.all_modules) >=1:
            self.on_rely(self.all_modules, self.set_output)
        else:
            self.all_modules[0].on("end", self.set_output)
        for module in self.all_modules:
            module.run()

    def set_output(self):
        if not os.path.exists(self.work_dir + "/tmp"):
            os.mkdir(self.work_dir + "/tmp")
        for module in self.all_modules:
            if len(os.listdir(module.output_dir)) >=1:
                link_dir(module.output_dir, self.work_dir + "/tmp")
        add_merge(self.work_dir + "/tmp", "result.xls", self.output_dir +"/all.coregene.xls")
        add_fasta(self.work_dir + "/tmp", "cor_gene.fa", self.output_dir +"/all.cor_gene.fa")
        self.end()

    def get_sample(self,file):
        list = []
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list.append(lin[0])
        return list

    def run(self):
        super(GetCoregeneModule, self).run()
        self.run_all_coregene()

    def end(self):
        super(GetCoregeneModule, self).end()