# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# version 1.0
# last_modify: 2021.04.28

import os,re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class ProdigalModule(Module):
    """
    小工具prodigal预测模块
    """

    def __init__(self, work_id):
        super(ProdigalModule, self).__init__(work_id)
        option = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的文件
            {"name": "gene_prefix", "type": "string", "default": "ORF"},  # 基因前缀，自定义前缀,默认gene
            {"name": "sample", "type": "string"},
            {"name": "cds_fnn", "type": "outfile", "format": "sequence.fasta"},  # 预测出的染色体序列
            {"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            {"name": "length_distribu", "type": "outfile", "format": "paternity_test.tab"},  # 样品基因序列的长度分布文件
            {"name": "sample_gene_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 样品编码基因预测gff文件
            {"name": "sample_trna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "sample_rrna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "sample_repeat_gff", "type": "outfile", "format": "gene_structure.gff3"},  # TRE预测结果
            {"name": "sample_gene_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample_gene_faa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "software_list", "type":"string", "default":"glimmer"},  #使用的软件，逗号分割
            {"name": "trans_code", "type" : "string","default":"11"},
            {"name": "p_trans_code","type": "string","default":"11"}, #
            {"name": "p_software_list", "type": "string","default":"genemark"},  #质粒的基因预测使用的软件，逗号分割


        ]
        self.add_option(option)
        self.cds_predict = self.add_module('tool_lab.gene_predicts')
        self.dna_predict_tidy = self.add_tool('tool_lab.dna_predict_tidy')

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genome").is_set:
            raise OptionError("必须提供基因组文件")
        if not self.option("sample"):
            raise OptionError("必须提供样品名")
        return True

    def run_cds_predict(self):
        opts = {
            'input_genome': self.option('genome'),
            'software_list': "prodigal",
            'orf_prefix': self.option('gene_prefix'),
            }
        opts['trans_code'] = self.option('trans_code') #zouguanqing
        self.cds_predict.set_options(opts)
        self.cds_predict.on('end', self.run_dna_predict_tidy)
        self.cds_predict.run()

    def run_dna_predict_tidy(self):
        if os.path.exists(self.cds_predict.output_dir + "/" + self.option("sample") + ".predict.gff"):
            with open(self.cds_predict.output_dir + "/" + self.option("sample") + ".predict.gff") as f:
                data = f.readlines()
                all_len = len(data)
            if all_len > 1:
                opts = {
                    'input_genome': self.option('genome'),
                    'gene_prefix': self.option('gene_prefix'),
                    'gene_gff': self.cds_predict.option('gff'),
                }
                self.dna_predict_tidy.set_options(opts)
                self.dna_predict_tidy.on('end', self.get_info)
                self.dna_predict_tidy.run()
        else:
            self.end()


    def get_info(self):
        if os.path.exists(self.output_dir + '/' + self.option("sample")):
            shutil.rmtree(self.output_dir + '/' + self.option("sample"))
        os.mkdir(self.output_dir + '/' + self.option("sample"))
        if os.path.exists(self.cds_predict.output_dir + "/" + self.option("sample") + ".predict.gff"):
            os.link(self.cds_predict.output_dir + "/" + self.option("sample") + ".predict.gff",self.output_dir + '/' + self.option("sample") + "/Prodigal_Gene prediction.gff3")
        if os.path.exists(self.cds_predict.output_dir + "/" + self.option("sample") + ".faa"):
            os.link(self.cds_predict.output_dir + "/" + self.option("sample") + ".faa", self.output_dir + '/' + self.option("sample") + "/Prodigal_Gene prediction.faa")
        if os.path.exists(self.cds_predict.output_dir + "/" + self.option("sample") + ".fnn"):
            os.link(self.cds_predict.output_dir + "/" + self.option("sample") + ".fnn", self.output_dir + '/' + self.option("sample") + "/Prodigal_Gene prediction.fna")
        for file in os.listdir(self.dna_predict_tidy.output_dir):
            if file.endswith(".gff"):
                with open(self.dna_predict_tidy.output_dir + "/" + file,"r") as f, open(self.output_dir + "/" + self.option("sample") + "/Prodigal_Gene prediction_detail.xls","w") as t:
                    t.write("Gene ID\tSample Name\tLocation\tStrand\tStart\tEnd\tGene Len（bp）\tInitiator Codon\tTerminator Codon\n")
                    data = f.readlines()
                    gene_num =0
                    gene_len = 0
                    for i in data[1:]:
                        tmp = i.strip().split("\t")
                        gene_num += 1
                        gene_len += abs(int(tmp[3]) - int(tmp[2]))
                        t.write(tmp[0] + "\t" + self.option("sample") + "\t" + tmp[1] + "\t" + tmp[4] + "\t" + tmp[2] +
                                "\t" + tmp[3] + "\t" + str(abs(int(tmp[3]) - int(tmp[2]))) + "\t" + tmp[9] + "\t" + tmp[10] + "\n")
        with open(self.option("genome").prop['path']) as v:
            GC_len = 0; N_len = 0; AT_len = 0
            for line in v:
                if line.startswith(">"):
                    pass
                else:
                    value = line.strip()
                    GC_len += (value.count('g') + value.count('G') + value.count('c') + value.count('C'))
                    N_len += (value.count('n') + value.count('N'))
                    AT_len += (value.count('a') + value.count('A') + value.count('t') + value.count('T'))
            all_len = GC_len + N_len + AT_len
            gc_percent = round(GC_len / float(all_len) *100,2)
        with open(self.output_dir + '/' + self.option("sample") + "/Prodigal_Gene prediction_stat.xls","w") as tt:
            tt.write("Sample Name\tGene No.\tGene Total Len（bp）\tGene Average Len（bp）\tGC Content  (%)\tGene/Genome (%)\n")
            if gene_num ==0:
                tt.write(self.option("sample") + "\t" + "0" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\n")
            else:
                tt.write(self.option("sample") + "\t" + str(gene_num) + "\t" + str(gene_len) + "\t" + str(round(float(gene_len)/gene_num *100,2)) +
                                "\t" + str(gc_percent) + "\t" + str(round(float(gene_len)/all_len *100,2)) + "\n")
        self.end()

    def run(self):
        super(ProdigalModule, self).run()
        self.run_cds_predict()

    def end(self):
        super(ProdigalModule, self).end()