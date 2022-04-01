# -*- coding: utf-8 -*-
# __author__ = "qiuping"
# last_modify_date:2017.02.14

from Bio import SeqIO
import numpy as np
from collections import defaultdict
from collections import Counter


class gene_dict(object):
    def __init__(self):
        self.len = 0
        self.tran = []
        self.full_name = ''
        self.seq = ''
        self.name = ''


class tran_dict(object):
    def __init__(self):
        self.len = 0
        self.name = ''
        self.seq = ''


def get_trinity_info(fasta):
    trinity_info = {}
    for seq in SeqIO.parse(fasta, "fasta"):
        tran = tran_dict()
        tran.name = seq.id
        tran.len = len(seq.seq)
        tran.seq = seq.seq
        if seq.id.split('_i')[0] not in trinity_info:
            gene = gene_dict()
            gene.name = seq.id.split('_i')[0]
            gene.full_name = seq.id
            gene.len = len(seq.seq)
            gene.seq = seq.seq
            gene.tran.append(tran)
            trinity_info[gene.name] = gene
        elif trinity_info[seq.id.split('_i')[0]].len >= tran.len:
            gene.tran.append(tran)
        else:
            gene.name = tran.name.split('_i')[0]
            gene.len = tran.len
            gene.seq = tran.seq
            gene.full_name = tran.name
            gene.tran.append(tran)

    return trinity_info


def stat_info(trinity_info, gene_path, stat_path, len_dir_path, length, full_name_path):
    gene_seq_num = len(trinity_info)
    tran_seq_num = 0
    gene_GC_num = 0
    tran_GC_num = 0
    gene_base_num = 0
    tran_base_num = 0
    gene_len = []
    tran_len = []
    gene_len_distri = {}
    tran_len_distri = {}
    for one in length.split(','):
        one = int(one)
        gene_len_distri[one] = defaultdict(int)
        tran_len_distri[one] = defaultdict(int)
    with open(gene_path, 'wb') as g, open(stat_path, 'wb') as s, open(full_name_path, 'wb') as fg:
        for i in trinity_info:
            gene = trinity_info[i]
            g.write('>{}\n{}\n'.format(i, gene.seq))
            fg.write('{}\n'.format(gene.full_name))
            tran_seq_num += len(gene.tran)
            gene_GC_num += (gene.seq.count('C') + gene.seq.count('G'))
            gene_base_num += gene.len
            gene_len.append(gene.len)
            for l in gene_len_distri:
                find_range(gene.len, l, gene_len_distri[l])
            for t in gene.tran:
                tran_GC_num += (t.seq.count('C') + t.seq.count('G'))
                tran_base_num += t.len
                tran_len.append(t.len)
                for l in tran_len_distri:
                    find_range(t.len, l, tran_len_distri[l])
        # gene_max_len, gene_min_len, gene_average_len, gene_N50, gene_N90 = len_stat(gene_len, gene_base_num)
        gene_len_stat = len_stat(gene_len, gene_base_num)
        # tran_max_len, tran_min_len, tran_average_len, tran_N50, tran_N90 = len_stat(tran_len, tran_base_num)
        tran_len_stat = len_stat(tran_len, tran_base_num)
        gene_GC_per = float(gene_GC_num) / gene_base_num * 100
        tran_GC_per = float(tran_GC_num) / tran_base_num * 100
        s.write('\tgenes\ttranscripts\tgene_count\ttranscript_count\ntotal seq num\t%s\t%s\t--\t--\ntotal base num\t%s\t%s\t--\t--\npercent GC\t%s\t%s\t--\t--\nlargest transcript\t%s\t%s\t%s\t%s\nsmallest transcript\t%s\t%s\t%s\t%s\naverage length\t%s\t%s\t%s\t%s\nN50\t%s\t%s\t%s\t%s\nN90\t%s\t%s\t%s\t%s' % (gene_seq_num, tran_seq_num, gene_base_num, tran_base_num, '%0.4g' % gene_GC_per, '%0.4g' % tran_GC_per, gene_len_stat['max'][0], tran_len_stat['max'][0], gene_len_stat['max'][1], tran_len_stat['max'][1], gene_len_stat['min'][0], tran_len_stat['min'][0], gene_len_stat['min'][1], tran_len_stat['min'][1], gene_len_stat['average'][0], tran_len_stat['average'][0], gene_len_stat['average'][1], tran_len_stat['average'][1], gene_len_stat['N50'][0], tran_len_stat['N50'][0], gene_len_stat['N50'][1], tran_len_stat['N50'][1], gene_len_stat['N90'][0], tran_len_stat['N90'][0], gene_len_stat['N90'][1], tran_len_stat['N90'][1]))
        for one in gene_len_distri:
            one = int(one)
            with open("{}/{}_length.distribut.txt".format(len_dir_path, one), "wb") as l:
                l.write('length\tgene_num\tgene_per\ttranscript_num\ttranscript_per\n')
                start = 1
                gene_steps = sorted(gene_len_distri[one].keys())
                tran_steps = sorted(tran_len_distri[one].keys())
                for i in range(one, tran_steps[-1] + 1, one):
                    if i in tran_steps:
                        tnum = tran_len_distri[one][i]
                        tper = '%0.4g' % (tran_len_distri[one][i] * 100.0 / tran_seq_num)
                    else:
                        tnum = 0
                        tper = 0
                    if i in gene_steps:
                        gnum = gene_len_distri[one][i]
                        gper = '%0.4g' % (gene_len_distri[one][i] * 100.0 / gene_seq_num)
                    else:
                        gnum = 0
                        gper = 0
                        # l.write('{}-{}\t{}\t{}\t{}\t{}\n'.format(start, i, gene_len_distri[one][i], '%0.4g' % (gene_len_distri[one][i] * 100.0 / gene_seq_num), tran_len_distri[one][i], '%0.4g' % (tran_len_distri[one][i] * 100.0 / tran_seq_num)))
                    l.write('{}-{}\t{}\t{}\t{}\t{}\n'.format(start, i, gnum, gper, tnum, tper))
                    start = i + 1
    return True


def find_range(len_, step, dict_):
    """
    计算某一个长度序列应该属于哪个区间，并将相应的dict 加1
    例如某条序列 长度len_为32，要计算步长20时，属于哪个区间，则传入参数应当是(32, 20, step_20)
    最后计算可知32 属于21-40的区间，字典step_20[40]应当加1

    :param len_:  序列的长度
    :param step:  步长
    :param dict_: 需要处理的字典
    """
    i = (len_ + step - 1) / step
    dict_[i * step] += 1


def len_stat(len_list, base_num):
    max_len = max(len_list)
    min_len = min(len_list)
    average_len = int(round(np.mean(len_list)))
    len_list.sort(reverse=True)
    len50 = 0
    len90 = 0
    N50 = 0
    N90 = 0
    for i in len_list:
        if len50 < base_num * 0.5:
            len50 += i
            N50 = i
        elif len90 < base_num * 0.9:
            len90 += i
            N90 = i
    len_count = Counter(len_list)
    stat_list = dict()
    stat_list.update(
        {
            'max': [max_len, len_count[max_len]],
            'min': [min_len, len_count[min_len]],
            'average': [average_len, len_count[average_len]],
            'N50': [N50, len_count[N50]],
            'N90': [N50, len_count[N90]]
        }
    )
    return stat_list

# a = get_trinity_info('C:\Users\ping.qiu.MAJORBIO\Desktop\output\\Trinity.fasta')
# stat_info(a, 'C:\Users\ping.qiu.MAJORBIO\Desktop\output\\gene.fasta', 'C:\Users\ping.qiu.MAJORBIO\Desktop\output\\trinity.fasta.stat.xls', 'C:\Users\ping.qiu.MAJORBIO\Desktop\output', '100,400', 'C:\Users\ping.qiu.MAJORBIO\Desktop\output\\gene_full_name.txt')
