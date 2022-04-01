# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import re
import os
import sys
import math
import pandas as pd


class Enrich(object):
    def __init__(self):
        self.diff_fc = {}
        self.max_gene_num = 100
        self.max_term_num = 1000
        self.min_go_dep = 3
        self.max_go_dep = 7
        self.max_p_value = 1
        self.max_padj_value = 1
        self.anno_list = []
        self.seq_id2name = dict()

    def get_fc(self, fc_table):
        '''
        读取表达量上下调信息
        '''
        fc_data = pd.read_table(fc_table, header=0)
        self.diff_fc = dict(zip(fc_data.seq_id, fc_data.log2fc))

    def set_gene_detail(self, gene_detail_file):
        '''
        根据基因描述信息获取基因或转录本与 gene name 对应关系
        '''
        with open(gene_detail_file, 'rb') as f:
            for line in f.readlines()[1:]:
                if line.split("\t")[2] != "":
                    self.seq_id2name[line.split("\t")[0]] = line.split("\t")[2]
                    self.seq_id2name[line.split("\t")[1]] = line.split("\t")[2]

    def filter_go_enrich(self, go_enrich_table):
        '''
        过滤go富集分析结果, 并得到富集弦图
        '''
        enrich_data = pd.read_table(go_enrich_table, header=0)

        # 筛选表格
        if self.anno_list:
            enrich_data = enrich_data[enrich_data.go_id.isin(self.anno_list)]
        else:
            enrich_data = enrich_data[enrich_data.p_uncorrected <= self.max_p_value]
            enrich_data = enrich_data[enrich_data.p_corrected <= self.max_padj_value]
            enrich_data = enrich_data[(enrich_data.depth >= self.min_go_dep) & (enrich_data.depth <= self.max_go_dep)]

        # 获取出现频率前 **的基因
        seqs = []
        seq_lists = list(enrich_data.seq_list)
        for i in seq_lists:
            seqs.extend(i.split("|"))
        counts = [[x, seqs.count(x)] for x in list(set(seqs))]
        counts = sorted(counts, key=lambda x: x[1], reverse=True)[:self.max_gene_num]
        choose_genes = [x[0] for x in counts]

        # 筛选基因
        choose_seq_lists_genes = [set(x.split("|")).intersection(set(choose_genes)) for x in seq_lists]
        enrich_data['choose_gene'] = pd.Series(choose_seq_lists_genes, index=enrich_data.index)
        choose_seq_lists_genes_num = [len(set(x.split("|")).intersection(set(choose_genes))) for x in seq_lists]
        enrich_data['choose_gene_num'] = pd.Series(choose_seq_lists_genes_num, index=enrich_data.index)
        enrich_data = enrich_data[(enrich_data.choose_gene_num > 0)]

        enrich_data = enrich_data[:self.max_term_num]

        choose_seqs = []
        choose_seq_lists = list(enrich_data.choose_gene)
        for i in choose_seq_lists:
            choose_seqs.extend(list(i))
        choose_seqs = list(set(choose_seqs))

        with open("go_enrich_choose.table", 'w') as f, open("go_enrich_detail.table", 'w') as f2:
            ids = list(enrich_data.go_id)
            terms = list(enrich_data.discription)
            f.write("{}\t{}\t{}\t{}\n".format("seq_id", "\t".join(ids), "log2FC", "Gene_name"))
            f2.write("{}\t{}\t{}\t{}\t{}\n".format("Seq", "GO_id",  "GO_description", "log2FC", "Gene_name"))
            for seq in choose_seqs:
                gene_name = seq
                if self.seq_id2name.has_key(seq):
                    gene_name = self.seq_id2name[seq]
                has_gene = []
                for line in range(0, len(enrich_data)):
                    annot_id = ids[line]
                    des = terms[line]
                    genes = list(list(enrich_data.choose_gene)[line])
                    if seq in genes and seq in self.diff_fc:
                        has_gene.append('1')
                        f2.write("{}\t{}\t{}\t{}\t{}\n".format(seq, annot_id, des, self.diff_fc[seq], gene_name))
                    else:
                        has_gene.append('0')
                f.write("{}\t{}\t{}\t{}\n".format(seq, "\t".join(has_gene), self.diff_fc[seq], gene_name))

        with open("enrich_zscore", 'w') as f:
            terms = list(enrich_data.discription)
            ids = list(enrich_data.go_id)
            seq_lists = list(enrich_data.seq_list)
            for i in range(0, len(terms)):
                fc_list = [self.diff_fc[x] for x in list(seq_lists[i].split("|")) if x in self.diff_fc ]
                fc_up = [self.diff_fc[x] for x in list(seq_lists[i].split("|")) if x in self.diff_fc and self.diff_fc[x] > 0]
                fc_down = [self.diff_fc[x] for x in list(seq_lists[i].split("|")) if x in self.diff_fc and self.diff_fc[x] < 0]
                z_score = (len(fc_up) - len(fc_down)) / (math.sqrt(len(fc_list))+0.01)
                f.write("{}\t{}\t{}\n".format(ids[i], terms[i], z_score))

    def filter_kegg_enrich(self, kegg_enrich_table):
        '''
        过滤kegg富集分析结果, 并得到富集弦图
        '''
        enrich_data = pd.read_table(kegg_enrich_table, header=0)
        if self.anno_list:
            enrich_data = enrich_data[enrich_data.id.isin(self.anno_list)]
        else:
            enrich_data = enrich_data[enrich_data.pvalue < self.max_p_value]
            enrich_data = enrich_data[enrich_data.corrected_pvalue < self.max_padj_value]
        print("筛选参数 p {} padj {} term_num {} terms {}".format(self.max_p_value, self.max_padj_value, self.max_term_num, self.anno_list))
        print(len(enrich_data))

        if len(enrich_data) == 0:
            return "筛选条件所得的富集结果为空请重新设置参数"

        # 获取出现频率前 50的基因
        seqs = []
        seq_lists = list(enrich_data.seq_list)
        for i in seq_lists:
            seqs.extend(i.split("|"))
        counts = [[x, seqs.count(x)] for x in list(set(seqs))]
        counts = sorted(counts, key=lambda x: x[1], reverse=True)[:self.max_gene_num]
        choose_genes = [x[0] for x in counts]

        # 筛选基因
        choose_seq_lists_genes = [set(x.split("|")).intersection(set(choose_genes)) for x in seq_lists]
        enrich_data['choose_gene'] = pd.Series(choose_seq_lists_genes, index=enrich_data.index)
        choose_seq_lists_genes_num = [len(set(x.split("|")).intersection(set(choose_genes))) for x in seq_lists]
        enrich_data['choose_gene_num'] = pd.Series(choose_seq_lists_genes_num, index=enrich_data.index)
        enrich_data = enrich_data[(enrich_data.choose_gene_num > 0)]

        enrich_data = enrich_data[:self.max_term_num]

        choose_seqs = []
        choose_seq_lists = list(enrich_data.choose_gene)
        for i in choose_seq_lists:
            choose_seqs.extend(list(i))
        choose_seqs = list(set(choose_seqs))

        with open("kegg_enrich_choose.table", 'w') as f, open("kegg_enrich_detail.table", 'w') as f2:
            ids = list(enrich_data.id)
            terms = list(enrich_data.term)
            f.write("{}\t{}\t{}\t{}\n".format("seq_id", "\t".join(ids), "log2FC", "Gene_name"))
            f2.write("{}\t{}\t{}\t{}\t{}\n".format("Gene","Pathway_id", "Pathway_description", "log2FC", "Gene_name"))
            for seq in choose_seqs:
                gene_name = seq
                if self.seq_id2name.has_key(seq):
                    gene_name = self.seq_id2name[seq]

                has_gene = []
                for line in range(0, len(enrich_data)):
                    annot_id = ids[line]
                    des = terms[line]
                    genes = list(list(enrich_data.choose_gene)[line])
                    if seq in genes and seq in self.diff_fc:
                        has_gene.append('1')
                        f2.write("{}\t{}\t{}\t{}\t{}\n".format(seq, annot_id, des, self.diff_fc[seq], gene_name))
                    else:
                        has_gene.append('0')
                if seq in self.diff_fc:
                    f.write("{}\t{}\t{}\t{}\n".format(seq, "\t".join(has_gene), self.diff_fc[seq], gene_name))

        with open("enrich_zscore", 'w') as f:
            terms = list(enrich_data.term)
            ids = list(enrich_data.id)
            seq_lists = list(enrich_data.seq_list)
            for i in range(0, len(terms)):
                fc_list = [self.diff_fc[x] for x in list(seq_lists[i].split("|")) if x in self.diff_fc]
                fc_up = [self.diff_fc[x] for x in list(seq_lists[i].split("|")) if x in self.diff_fc and self.diff_fc[x] > 0]
                fc_down = [self.diff_fc[x] for x in list(seq_lists[i].split("|")) if x in self.diff_fc and self.diff_fc[x] < 0]
                z_score = (len(fc_up) - len(fc_down)) / (math.sqrt(len(fc_list))+0.01)
                f.write("{}\t{}\t{}\n".format(ids[i], terms[i], z_score))


if __name__ == "__main__":
    enrich = Enrich()
    enrich.get_fc(sys.argv[1])
    # enrich.filter_go_enrich(sys.argv[2])
    enrich.filter_kegg_enrich(sys.argv[2])
