#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,time
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module

class HousegeneModule(Module):
    """
    持家基因预测
    """

    def __init__(self, work_id):
        super(HousegeneModule, self).__init__(work_id)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"}
        ]
        self.add_option(options)
        self.predict_seq = self.add_module("bac_comp_genome.gene_predict_seq")
        self.core_gene = self.add_tool('bac_comp_genome.get_core_gene')
        #self.sample_name = self.option("fasta").prop['path'].split("/")[-1].split(".")[0]
        #self.modules = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option('fasta').is_set:
            raise OptionError("请传入序列文件！", code="")

    def run_gene_predict(self):
        opts = ({
            "genome": self.option("fasta"),
            "genome_name" : self.option("sample_name"),
            "gene_prefix" : "gene",
            "genome_type" : "draft"
        })
        self.predict_seq.set_options(opts)
        self.predict_seq.on("end", self.run_coregene)
        self.predict_seq.run()

    def run_coregene(self):
        opts = {
            "seq_faa": self.predict_seq.output_dir + "/"  + self.option("sample_name") + "_CDS.faa",
            "seq_gff": self.predict_seq.output_dir + "/" + self.option("sample_name") + "_CDS.gff",
            "sample_name": self.option("sample_name"),
            "method": "genome",
        }
        self.core_gene.set_options(opts)
        self.core_gene.on("end", self.set_output)
        self.core_gene.run()

    def run(self):
        super(HousegeneModule, self).run()
        self.run_gene_predict()

    def set_output(self):
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + "/" + self.option("sample_name")):
            shutil.rmtree(self.output_dir + "/" + self.option("sample_name"))
        os.mkdir(self.output_dir + "/" + self.option("sample_name"))
        os.link(self.core_gene.output_dir + "/" + self.option("sample_name") + '.cor_gene.fa', self.output_dir + "/" + self.option("sample_name") + "/" + self.option("sample_name") + "_all.hgene.fa")
        raw_hgene_file = self.core_gene.output_dir + "/" + self.option("sample_name") + ".result.xls"
        hgenefa_file = self.output_dir + "/" + self.option("sample_name") + "/" + self.option("sample_name") + "_hgene.fa"
        blast_file = self.output_dir + "/" + self.option("sample_name") + "/" + self.option("sample_name") + "_hgene_blast.xls"
        stat_file = self.output_dir + "/" + "stat_file.txt"
        blast_m8_file = self.core_gene.work_dir + "/" + self.option("sample_name") + ".matches.m8"
        with open(raw_hgene_file,"r") as r, open(hgenefa_file,"w") as h, open(blast_file,"w") as b, open(stat_file,"w") as s,open(blast_m8_file,"r") as blast:
            data1 = r.readlines()
            data2 = blast.readlines()
            s.write(self.option("sample_name") + "\t" + str(len(data1)-1))
            b.write("Genome" + "\t" + "Gene ID" + "\t" + "Location" + "\t" + "Start" + "\t" + "End" + "\t" + "Name" + "\t" + "Indentity" + "\t" + "Coverage" + "\t" + "Evalue" + "\t" + "Score""\n")
            for i in data1[1:]:
                for x in data2:
                    if  i.strip().split("\t")[5] == x.strip().split("\t")[0]:
                        b.write("\t".join(i.split("\t")[0:8]) + "\t" + x.strip().split("\t")[10] + "\t" + x.strip().split("\t")[11] + "\n")
                        h.write(">" + i.split("\t")[0] + "|" + i.split("\t")[1] + "|" + i.split("\t")[2] + ":" + i.split("\t")[3] + "-" + i.split("\t")[4] + "|" + i.split("\t")[5] + "\n" + i.split("\t")[8])
        if os.path.getsize(hgenefa_file) == 0:
            shutil.rmtree(self.output_dir + "/" + self.option("sample_name"))
        self.end()

    def end(self):
        super(HousegeneModule, self).end()