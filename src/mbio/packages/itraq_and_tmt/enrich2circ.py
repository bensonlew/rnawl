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
        self.max_gene_num = 50
        self.max_term_num = 15
        self.min_go_dep = 2
        self.max_go_dep = 15
        self.max_p_value = 1 

    def get_fc(self, fc_table):
        '''
        读取表达量上下调信息
        '''
        fc_data = pd.read_table(fc_table, header=0, dtype = {0:str})
        self.diff_fc = dict(zip(fc_data.accession_id,fc_data.log2fc))

    def filter_go_enrich(self, go_enrich_table):
        '''
        过滤go富集分析结果, 并得到富集弦图
        '''
        enrich_data = pd.read_table(go_enrich_table, header=0)
        enrich_data = enrich_data[enrich_data.p_uncorrected < self.max_p_value]
        enrich_data = enrich_data[(enrich_data.depth >= self.min_go_dep) & (enrich_data.depth <= self.max_go_dep)]
        
        # 获取出现频率前 50的基因
        seqs = []
        seq_lists=list(enrich_data.seq_list)
        for i in seq_lists:
            seqs.extend(i.split(";"))
        counts=[[x,seqs.count(x)] for x in list(set(seqs))]
        counts=sorted(counts, key=lambda x:x[1], reverse=True)[:self.max_gene_num]
        choose_genes = [x[0] for x in counts]
        
        # 筛选基因
        choose_seq_lists_genes = [set(x.split(";")).intersection(set(choose_genes)) for x in seq_lists]
        enrich_data['choose_gene'] = pd.Series(choose_seq_lists_genes, index=enrich_data.index)
        choose_seq_lists_genes_num = [len(set(x.split(";")).intersection(set(choose_genes))) for x in seq_lists]
        enrich_data['choose_gene_num'] = pd.Series(choose_seq_lists_genes_num, index=enrich_data.index)
        enrich_data = enrich_data[(enrich_data.choose_gene_num > 0)]
        
        enrich_data = enrich_data[:self.max_term_num]
        
        choose_seqs = []
        choose_seq_lists=list(enrich_data.choose_gene)
        for i in choose_seq_lists:
            choose_seqs.extend(list(i))
        choose_seqs = list(set(choose_seqs))
        
        with open("go_enrich_choose.table", 'w') as f , open("go_enrich_detail.table", 'w') as f2:
            ids =  list(enrich_data.go_id)
            terms = list(enrich_data.discription)
            f.write("{}\t{}\t{}\n".format("seq_id", "\t".join(ids), "log2FC"))
            f2.write("{}\t{}\t{}\n".format("Acession", "GO", "log2FC"))
            for seq in choose_seqs:
                has_gene = []
                for line in range(0, len(enrich_data)):
                    des = terms[line]
                    genes = list(list(enrich_data.choose_gene)[line])
                    if seq in genes:
                        has_gene.append('1')
                        f2.write("{}\t{}\t{}\n".format(seq, des, self.diff_fc[str(seq)]))
                    else:
                        has_gene.append('0')
                f.write("{}\t{}\t{}\n".format(seq, "\t".join(has_gene), self.diff_fc[str(seq)]))
        
        with open("enrich_zscore", 'w') as f:
            terms = list(enrich_data.discription)
            ids = list(enrich_data.go_id)
            seq_lists =  list(enrich_data.seq_list)
            for i in range(0, len(terms)):
                fc_list = [self.diff_fc[x] for x in list(seq_lists[i].split(";"))]
                fc_up = [self.diff_fc[x] for x in  list(seq_lists[i].split(";")) if self.diff_fc[x] > 0]
                fc_down = [self.diff_fc[x] for x in  list(seq_lists[i].split(";")) if self.diff_fc[x] < 0]
                z_score = (len(fc_up) - len(fc_down))/math.sqrt(len(fc_list))
                f.write("{}\t{}\t{}\n".format(ids[i], terms[i], z_score))

    def filter_kegg_enrich(self, kegg_enrich_table):
        '''
        过滤kegg富集分析结果, 并得到富集弦图
        '''
        enrich_data = pd.read_table(kegg_enrich_table, header=0)
        enrich_data = enrich_data[enrich_data.pvalue < self.max_p_value]
        
        # 获取出现频率前 50的基因
        seqs = []
        seq_lists=list(enrich_data.seq_list)
        for i in seq_lists:
            seqs.extend(i.split("|"))
        counts=[[x,seqs.count(x)] for x in list(set(seqs))]
        counts=sorted(counts, key=lambda x:x[1], reverse=True)[:self.max_gene_num]
        choose_genes = [x[0] for x in counts]
        
        # 筛选基因
        choose_seq_lists_genes = [set(x.split("|")).intersection(set(choose_genes)) for x in seq_lists]
        enrich_data['choose_gene'] = pd.Series(choose_seq_lists_genes, index=enrich_data.index)
        choose_seq_lists_genes_num = [len(set(x.split("|")).intersection(set(choose_genes))) for x in seq_lists]
        enrich_data['choose_gene_num'] = pd.Series(choose_seq_lists_genes_num, index=enrich_data.index)
        enrich_data = enrich_data[(enrich_data.choose_gene_num > 0)]
        
        enrich_data = enrich_data[:self.max_term_num]
        
        choose_seqs = []
        choose_seq_lists=list(enrich_data.choose_gene)
        for i in choose_seq_lists:
            choose_seqs.extend(list(i))
        choose_seqs = list(set(choose_seqs))
        
        with open("kegg_enrich_choose.table", 'w') as f , open("kegg_enrich_detail.table", 'w') as f2:
            ids =  list(enrich_data.id)
            terms = list(enrich_data.term)
            f.write("{}\t{}\t{}\n".format("seq_id", "\t".join(ids), "log2FC"))
            f2.write("{}\t{}\t{}\n".format("Acession", "KEGG", "log2FC"))
            for seq in choose_seqs:
                has_gene = []
                for line in range(0, len(enrich_data)):
                    des = terms[line]
                    genes = list(list(enrich_data.choose_gene)[line])
                    if seq in genes:
                        has_gene.append('1')
                        f2.write("{}\t{}\t{}\n".format(seq, des, self.diff_fc[seq]))
                    else:
                        has_gene.append('0')
                f.write("{}\t{}\t{}\n".format(seq, "\t".join(has_gene), self.diff_fc[seq]))
        
        with open("enrich_zscore", 'w') as f:
            terms = list(enrich_data.term)
            ids = list(enrich_data.id)
            seq_lists =  list(enrich_data.seq_list)
            for i in range(0, len(terms)):
                fc_list = [self.diff_fc[x] for x in list(seq_lists[i].split("|"))]
                fc_up = [self.diff_fc[x] for x in  list(seq_lists[i].split("|")) if self.diff_fc[x] > 0]
                fc_down = [self.diff_fc[x] for x in  list(seq_lists[i].split("|")) if self.diff_fc[x] < 0]
                z_score = (len(fc_up) - len(fc_down))/math.sqrt(len(fc_list))
                f.write("{}\t{}\t{}\n".format(ids[i], terms[i], z_score))


if __name__ == "__main__":
    enrich = Enrich()
    enrich.get_fc(sys.argv[1])
    # enrich.filter_go_enrich(sys.argv[2])
    enrich.filter_kegg_enrich(sys.argv[2])
