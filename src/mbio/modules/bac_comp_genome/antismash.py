#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,re
import shutil
import gevent
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class AntismashModule(Module):
    """
    单个基因组次级代谢产物的预测
    last_modify: 2019.10.10
    """
    def __init__(self, work_id):
        super(AntismashModule, self).__init__(work_id)
        options = [
            {"name": "gbk_dir", "type": "infile", "format": "gene_structure.gbk_dir"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},
        ]
        self.step.add_steps('antismash')
        self.get_fasta = self.add_tool("bac_comp_genome.antismash_fasta")
        self.add_option(options)
        self.modules = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option('gbk_dir').is_set:
            raise OptionError("请设置基因组基因gbk文件夹不存在！")
        if not self.option('sample_name'):
            raise OptionError("请设置样品名称！")
        if not self.option('genome'):
            raise OptionError("请设置基因组序列文件！")

    def run_antismash(self):
        """
        antismash 运行
        :return:
        """
        if self.option('gbk_dir').is_set:
            files = os.listdir(self.option('gbk_dir').prop['path'])
            for file in files:
                if re.search(r'.gbk', file):
                    n = self.get_size(self.option('gbk_dir').prop['path'] + '/' + file)
                    if int(n) >= 10000:
                        self.antismash = self.add_tool('bac_comp_genome.antismash')
                        opts = {
                            "genome_gbk": self.option('gbk_dir').prop['path'] + '/' + file,
                        }
                        self.antismash.set_options(opts)
                        self.modules.append(self.antismash)
            if len(self.modules) > 1:
                self.on_rely(self.modules, self.run_get_fasta)
            elif len(self.modules) == 1:
                self.modules[0].on('end', self.run_get_fasta)
            for module in self.modules:
                module.run()

    def get_size(self, gbk):
        with open (gbk, "r") as f:
            line =f.readlines()[0]
            m = re.search("(\s+[0-9]*)\sbp", line)
            num = m.group(1)
        return num

    def get_summary(self):
        if os.path.exists(self.output_dir + "/" + self.option("sample_name")):
            shutil.rmtree(self.output_dir + "/" + self.option("sample_name"))
        os.mkdir(self.output_dir + "/" + self.option("sample_name"))
        n=0
        with open (self.work_dir + "/" + self.option('sample_name') + ".antismash_anno.xls", "w") as g:
            g.write("Cluster ID\tLocation\tType\tStart\tEnd	\tGene No.\tGenes\tpredicted_structure\tMost Similar Cluster\tSimilarity\n")
            for module in self.modules:
                if os.path.exists(module.output_dir + "/antismash_anno.xls"):
                    with open (module.output_dir + "/antismash_anno.xls", "r") as f:
                        lines = f.readlines()
                        for line in lines[1:]:
                            n +=1
                            lin = line.strip().split("\t")
                            if lin[7] != "-":
                                self.logger.info("aaa")
                                self.logger.info(module.output_dir + "/core_structures/"+ lin[0].capitalize() + ".png")
                                os.link(module.output_dir + "/core_structures/" + lin[0].lower() + ".png", self.output_dir + "/" + self.option("sample_name") + "/" + "Cluster"+ str(n) + ".png")
                            des = "Cluster"+ str(n) + "\t" + lin[10] + "\t" +"\t".join(lin[1:10])
                            g.write(des + "\n")

    def run_get_fasta(self):
        self.logger.info("sssssssssssss1")
        self.get_summary()
        self.logger.info("sssssssssssss2")
        self.num =self.get_num(self.work_dir + "/" + self.option('sample_name') + ".antismash_anno.xls")
        if self.num >1:
            self.get_fasta.set_options({
                "sample_name": self.option("sample_name"),
                "genome": self.option("genome"),
                "antismash": self.work_dir + "/" + self.option('sample_name') + ".antismash_anno.xls"
            })
            self.get_fasta.on("end", self.set_output)
            self.get_fasta.run()
        else:
            shutil.rmtree(self.output_dir + "/" + self.option("sample_name"))
            self.end()

    def run(self):
        super(AntismashModule, self).run()
        self.run_antismash()

    def get_num(self, file):
        with open(file, "r") as f:
            lines = f.readlines()
        return len(lines)

    def set_output(self):
        link_file(self.work_dir + "/" + self.option('sample_name') + ".antismash_anno.xls", self.output_dir + "/"  + self.option('sample_name') + ".antismash_anno.xls")
        if not len(os.listdir(self.output_dir + "/" + self.option("sample_name"))) > 0:
            shutil.rmtree(self.output_dir + "/" + self.option("sample_name"))
        link_dir(self.get_fasta.output_dir, self.output_dir)
        self.end()

    def end(self):
        super(AntismashModule, self).end()