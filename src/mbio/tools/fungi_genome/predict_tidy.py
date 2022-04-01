# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# version 1.0
# last_modify: 2018.02.28

import os
import re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import Counter


class PredictTidyAgent(Agent):
    """
    对基因预测、tRNA预测和rRNA预测的结果进行整理：
    1. 利用tRNA结果去除假基因（ORF预测基因中与tRNA预测结果重叠的部分）；
    2. 将预测结果按位置信息进行统一排序；
    3. 统一基因id：默认为gene0001开始，可设置前缀替换“gene”
    """

    def __init__(self, parent):
        super(PredictTidyAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测结果
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "gene_prefix", "type": "string", "default": "gene"},  # 基因前缀，自定义前缀,默认gene
            {"name": "plasmid_prefix", "type": "string", "default": ""},
            # 质粒前缀,字典形式含一个或多个质粒，"{‘plassmidA’: ‘plasa', ‘plassmidB’: ‘plasb'}"
            {"name": "gene", "type": "outfile", "format": "gene_structure.gff3"},  # 基因预测结果,排序整理后的，并加上累积长度起始位置
            {"name": "trna", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "rrna", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "ffn", "type": "string", "default": ""},
            {"name": "faa", "type": "string", "default": ""},
            {"name": "trna_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "rrna_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample", "type": "string", "default": "sample"}

        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gene_gff").is_set:
            raise OptionError("必须设置参数gene_gff", code="32101801")
        if not self.option("trna_gff").is_set:
            raise OptionError("必须设置参数trna_gff", code="32101802")
        if not self.option("rrna_gff").is_set:
            raise OptionError("必须设置参数rrna_gff", code="32101803")
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome", code="32101804")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(PredictTidyAgent, self).end()


class PredictTidyTool(Tool):
    def __init__(self, config):
        super(PredictTidyTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]

        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"

    def sort_number(self):
        gene_table_path = self.option("gene_gff").prop["path"]
        gene_table = pd.read_table(gene_table_path, sep='\t', header=0)
        gene = gene_table.ix[:, ["Sequence id", "Start", "End"]]
        if len(gene_table) > 0:
            gene["location"] = gene["Sequence id"].str.split("_", expand=True)[0]
            gene_sort = gene["location"].str.split("caffold", expand=True)
        gene["s"] = gene.ix[:, ["Start", "End"]].apply(lambda x: x.min(), axis=1)



        trna_table_path = self.option("trna_gff").prop["path"]
        trna_table = pd.read_table(trna_table_path, sep='\t', header=0)
        trna = trna_table.ix[:, ["Sequence id", "Start", "End"]]
        trna["s"] = trna.ix[:, ["Start", "End"]].apply(lambda x: x.min(), axis=1)
        #trna["e"] = trna.ix[:, ["Start", "End"]].apply(lambda x: x.max(), axis=1)
        if len(trna_table) > 0:
            trna["location"] = trna["Sequence id"].str.split("_", expand=True)[0]
            trna_sort = trna["location"].str.split("caffold", expand=True)
        rrna_table_path = self.option("rrna_gff").prop["path"]
        rrna_table = pd.read_table(rrna_table_path, sep='\t', header=0)
        rrna = rrna_table.ix[:, ["Sequence id", "Start", "End"]]
        rrna["s"] = rrna.ix[:, ["Start", "End"]].apply(lambda x: x.min(), axis=1)
        if len(rrna_table) > 0:
            rrna["location"] = rrna["Sequence id"].str.split("_", expand=True)[0]
            rrna_sort = rrna["location"].str.split("caffold", expand=True)
        if re.match("^[A-Za-z]caffold[0-9]+$", list(gene["location"])[1]):
            if len(trna_table) > 0:
                trna["sort"] = trna_sort[1].astype(int)
            if len(gene_table) > 0:
                gene["sort"] = gene_sort[1].astype(int)
            if len(rrna_table) > 0:
                rrna["sort"] = rrna_sort[1].astype(int)
        else:
            if len(trna_table) > 0:
                trna["sort"] = trna_sort[0]
            if len(gene_table) > 0:
                gene["sort"] = gene_sort[0]
            if len(rrna_table) > 0:
                rrna["sort"] = rrna_sort[0]
        summary = pd.concat([gene.ix[:, ["Sequence id", "s", "sort"]], trna.ix[:, ["Sequence id", "s", "sort"]],
                             rrna.ix[:, ["Sequence id", "s", "sort"]]])
        summary = summary.sort_values(by=["sort", "s"])

        gene_id_list = []
        if self.option("plasmid_prefix") != "":
            pre = eval(self.option("plasmid_prefix"))  # 如果字符串引号为双引号，则通不过该不
            for k in sorted(pre.keys()):
                plas_one = summary[summary["sort"] == k]
                gene_suffix_len = len(str(len(plas_one)))
                if gene_suffix_len < 5:
                    gene_suffix_len = 4
                for n in range(len(plas_one)):
                    p = n + 1
                    gene_id = pre[k] + "0" * (gene_suffix_len - len(str(p))) + str(p)
                    gene_id_list.append(gene_id)
        else:
            gene_suffix_len = len(str(len(summary)))
            if gene_suffix_len < 5:
                gene_suffix_len = 4
            for n in range(len(summary)):
                k = n + 1
                gene_id = self.option("gene_prefix") + "0" * (gene_suffix_len - len(str(k))) + str(k)
                gene_id_list.append(gene_id)
        summary["Gene ID"] = gene_id_list
        del summary["s"]
        del summary["sort"]
        # "提取排序好的gene_id, 形成各预测的gff文件"
        del gene_table["Gene id"]
        gene = summary.merge(gene_table, on='Sequence id', how='inner')

        self.seq_gene_dic = {}
        self.gene_sort_list = []
        for i in range(len(gene)):
            self.seq_gene_dic[gene["Sequence id"][i]] = gene["Gene ID"][i]
            self.gene_sort_list.append(gene["Sequence id"][i])
        #del gene["Gene id"]
        gene_file_path = self.work_dir + "/predict.gff"
        gene = gene.set_index("Gene ID")
        # if len(gene) > 0:
        #     location = gene["Sequence id"].str.split("_", expand=True)[0]
        #     loc_nu = (Counter(location))
        #     for i in range(len(gene)):
        #         seq_id = gene["Sequence id"][i]
        #         gene_suffix_len = len(str(loc_nu[seq_id.split("_ORF")[0]]))
        #         nu = seq_id.split("_ORF")[1].lstrip("0")
        #         tmp = "_ORF" + "0" * (gene_suffix_len - len(nu)) + nu
        #         gene["Sequence id"][i] = re.subn("_ORF[0-9]*", tmp, gene["Sequence id"][i])[0]
        gene.to_csv(gene_file_path, sep="\t")



        rrna = summary.merge(rrna_table, on='Sequence id', how='inner')
        del rrna["Sequence Name"]
        del rrna["Type"]
        rrna_file_path = self.work_dir + "/rRNA.gff"
        rrna = rrna.set_index("Gene ID")
        if len(rrna) > 0:
            #location = rrna["Sequence id"].str.split("_", expand=True)[0]
            #loc_nu = (Counter(location))
            nu = 1
            gene_suffix_len = 3
            for i in range(len(rrna)):
                seq_id = rrna["Sequence id"][i]
                #gene_suffix_len = len(str(loc_nu[seq_id.split("_rRNA")[0]]))
                #nu = seq_id.split("_rRNA")[1].lstrip("0")
                tmp = "_rRNA" + "0" * (gene_suffix_len - len(str(nu))) + str(nu)
                rrna["Sequence id"][i] = re.subn("_rRNA[0-9]*", tmp, rrna["Sequence id"][i])[0]
                nu += 1
        rrna.to_csv(rrna_file_path, sep="\t")
        trna = summary.merge(trna_table, on='Sequence id', how='inner')
        del trna["Sequence Name"]
        trna_file_path = self.work_dir + "/tRNA.gff"
        trna = trna.set_index("Gene ID")
        if len(trna) > 0:

            nu = 1
            gene_suffix_len = 3
            for i in range(len(trna)):
                seq_id = trna["Sequence id"][i]
                #gene_suffix_len = len(str(loc_nu[seq_id.split("_tRNA")[0]]))
                #nu = seq_id.split("_tRNA")[1].lstrip("0")

                tmp = "_tRNA" + "0" * (gene_suffix_len - len(str(nu))) + str(nu)
                trna["Sequence id"][i] = re.subn("_tRNA[0-9]*", tmp, trna["Sequence id"][i])[0]
                nu += 1
        trna.to_csv(trna_file_path, sep="\t")


    def run_get_fa_bygff(self, genome, gff):
        cmd = '{} {}get_fa_bygff.pl {} {} {} {}'.format(self.perl_path, self.perl_script, genome,
                                                        gff, "0", self.work_dir)
        name = gff.split('/')[-1].split(".")[0]
        command = self.add_command(str(name.lower()), cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn和整理GFF文件运行完成")
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!", code="32101801")


    def fasta_change_id(self, fasta, seq_gene_dic, outname):
        fr = open(fasta)
        fw = open(self.output_dir+'/'+outname, 'w')
        fa_dic = {}
        for line in fr:
            if line.startswith('>'):
                sp_line = line.split(' ')
                if len(sp_line) < 4:
                    #fw.write(line)
                    k = line
                    fa_dic[k] = line
                    continue
                k = sp_line[3]
                if k in seq_gene_dic.keys():
                    fa_dic[k]='>'+seq_gene_dic[k]+' '+' '.join(sp_line[1:])
                else:
                    fa_dic[k]= line
                    self.logger.info(k+" 在基因gff文件中不存在。该序列没有输出到fasta文件")
                #fw.write('>'+seq_gene_dic[k]+' '+' '.join(sp_line[1:]))
            else:
                fa_dic[k] += line
                #fw.write(line)

        for k in self.gene_sort_list:
            fw.write(fa_dic[k])



    def set_output(self):
        old_files = "{0}/predict.gff,{0}/tRNA.gff,{0}/rRNA.gff,{0}/{1}_rRNA.fnn,{0}/{1}_tRNA.fnn".format(self.work_dir, self.sample_name).split(',')
        new_files = "{0}/{1}_CDS.gff,{0}/{1}_tRNA.gff,{0}/{1}_rRNA.gff,{0}/{1}_rRNA.fnn,{0}/{1}_tRNA.fnn".format(self.output_dir, self.option('sample')).split(',')
        for i in new_files:
            if os.path.exists(i):
                os.remove(i)
            os.link(old_files[new_files.index(i)],i)

        #self.option('gene', self.output_dir + "/" + self.option('sample') + "_CDS.gff")
        #self.option('trna', self.output_dir + "/" + self.sample_name + "_tRNA.gff")
        #self.option('rrna', self.output_dir + "/" + self.sample_name + "_rRNA.gff")
        # if os.path.getsize(self.output_dir + "/" + self.sample_name + "_predict.fnn"):
        #     self.option('gene_fnn', self.output_dir + "/" + self.sample_name + "_predict.fnn")
        #     self.option('gene_faa', self.output_dir + "/" + self.sample_name + "_predict.faa")
        #if os.path.getsize(self.output_dir + "/" + self.sample_name + "_tRNA.fnn"):
        #    self.option('trna_fnn', self.output_dir + "/" + self.sample_name + "_tRNA.fnn")
        #if os.path.getsize(self.output_dir + "/" + self.sample_name + "_rRNA.fnn"):
        #    self.option('rrna_fnn', self.output_dir + "/" + self.sample_name + "_rRNA.fnn")

    def run(self):
        super(PredictTidyTool, self).run()
        self.sort_number()
        for i in ["/tRNA.gff", "/rRNA.gff"]:
            self.run_get_fa_bygff(self.option("input_genome").prop["path"], self.work_dir + i)
        self.fasta_change_id(self.option('ffn'), self.seq_gene_dic, self.option('sample') + "_CDS.fnn")
        self.fasta_change_id(self.option('faa'), self.seq_gene_dic, self.option('sample') + "_CDS.faa")

        self.set_output()
        self.end()
