# !/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "qiuping"
# last_modify:20160616

import math
import numpy as np


def read_matrix(matrixfile, control, other, replicates=None):
    """传入计数或表达量矩阵，返回一个字典，键为gene_id,值为对照组（样本）和实验组（样本）的count值，fpkm值；control为对照组名，other为实验组名"""
    with open(matrixfile, 'rb') as f:
        f_info = f.readlines()
        sample = f_info[0].strip('\n').split('\t')
        info_dict = {}
        if replicates == None:
            index_1 = sample.index(control)
            index_2 = sample.index(other)
            for line in f_info[1:]:
                aline = line.strip('\n').split('\t')
                info_dict[aline[0]] = []
                info_dict[aline[0]].append(aline[index_1])
                info_dict[aline[0]].append(aline[index_2])
            return info_dict
        else:
            with open(replicates, 'rb') as g:
                group_dict = dict()
                group_dict[control] = []
                group_dict[other] = []
                for line in g.readlines():
                    line = line.strip('\n').split('\t')
                    if line[0] == control:
                        group_dict[control].append(line[1])
                    elif line[0] == other:
                        group_dict[other].append(line[1])
                group_sample = group_dict[control] + group_dict[other]
            for line in f_info[1:]:
                aline = line.strip('\n').split('\t')
                control_stat = []
                other_stat = []
                for sam in group_dict[control]:
                    control_stat.append(float(aline[sample.index(sam)]))
                for sam in group_dict[other]:
                    other_stat.append(float(aline[sample.index(sam)]))
                info_dict[aline[0]] = []
                info_dict[aline[0]].append(control_stat)
                info_dict[aline[0]].append(other_stat)
            return info_dict, group_sample


def mean_fpkm(count_dict, fpkm_dict):
    mean_fpkm = dict()
    control_libsize = 0
    other_libsize = 0
    con_count = []
    oth_count = []
    for gene in count_dict.keys():
        control_libsize += sum(count_dict[gene][0])
        other_libsize += sum(count_dict[gene][1])
        con_count.append(count_dict[gene][0])
        oth_count.append(count_dict[gene][1])
    c_counts = np.array(con_count).sum(axis=0)  # 计算对照组每个样本的count的总和
    o_counts = np.array(oth_count).sum(axis=0)  # 计算实验组每个样本的count的总和
    for gene in fpkm_dict.keys():
        mean_fpkm[gene] = []
        con = 0
        oth = 0
        for i in range(len(fpkm_dict[gene][0])):
            con += fpkm_dict[gene][0][i] * c_counts[i]
            oth += fpkm_dict[gene][1][i] * o_counts[i]
        c_mean = '%0.3f' % (con / control_libsize)
        o_mean = '%0.3f' % (oth / other_libsize)
        mean_fpkm[gene] = [c_mean, o_mean]
    return mean_fpkm


def stat_edger(edgr_result, countfile, fpkmfile, control, other, output, replicates=None, diff_ci=0.05, regulate=True):
    """对edgeR结果进行统计，获得两两分组（或样本）的edgeR统计文件以及差异基因列表文件"""
    with open(edgr_result, 'rb') as e, open('%s/gene_count.txt.%s_vs_%s.edgr_stat.xls' % (output, control, other), 'wb') as w:
        eline = e.readline()
        if replicates is None:
            count_dict = read_matrix(countfile, control, other)
            fpkm_dict = read_matrix(fpkmfile, control, other)
            w.write("Gene_id\t%s_count\t%s_conut\t%s_fpkm\t%s_fpkm\tLog2FC(%s/%s)\tPvalue\tFDR\tSignificant\tRegulate\n" % (control, other, control, other, other, control))
        else:
            count_dict, group_sample = read_matrix(countfile, control, other, replicates)
            fpkm_dict, group_sample = read_matrix(fpkmfile, control, other, replicates)
            fpkm_mean = mean_fpkm(count_dict, fpkm_dict)
            head = "Gene_id\t"
            for sam in group_sample:
                head += "%s_count\t%s_fpkm\t" % (sam, sam)
            head += "%s_mean_fpkm\t%s_mean_fpkm\tLog2FC(%s/%s)\tPvalue\tFDR\tSignificant\tRegulate\n" % (control, other, other, control)
            w.write(head)
        # print '######%s-%s:regulate is %s' %  (control, other,regulate)
        while True:
            eline = e.readline().strip('\n').split('\t')
            if not eline[0]:
                break
            gene_id = eline[0]
            if replicates is None:
                fc = (float(fpkm_dict[gene_id][1]) + 0.1) / (float(fpkm_dict[gene_id][0]) + 0.1)
            else:
                fc = (float(fpkm_mean[gene_id][1]) + 0.1) / (float(fpkm_mean[gene_id][0]) + 0.1)
            logfc = math.log(fc, 2)
            if regulate:
                if logfc > 0:
                    reg = 'up'
                elif logfc < 0:
                    reg = 'down'
                else:
                    reg = 'no change'
            else:
                reg = 'undone'
            fdr = float(eline[4])
            if fdr < diff_ci:
                sig = 'yes'
            else:
                sig = 'no'
            if replicates is None:
                w.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_id, count_dict[gene_id][0], count_dict[gene_id][1], fpkm_dict[gene_id][0], fpkm_dict[gene_id][1], '%0.3f' % logfc, '%0.3f' % float(eline[3]), '%0.3f' % fdr, sig, reg))
                fpkm_dict.pop(gene_id)
            else:
                w_text = "%s\t" % gene_id
                count_list = count_dict[gene_id][0] + count_dict[gene_id][1]
                fpkm_list = fpkm_dict[gene_id][0] + fpkm_dict[gene_id][1]
                for i in range(len(count_list)):
                    w_text += "%s\t%s\t" % (count_list[i], fpkm_list[i])
                w.write("%s%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (w_text, fpkm_mean[gene_id][0], fpkm_mean[gene_id][1], '%0.3f' % logfc, '%0.3f' % float(eline[3]), '%0.3f' % fdr, sig, reg))
                fpkm_dict.pop(gene_id)
        if replicates == None:
            for gene in fpkm_dict.keys():
                w.write("%s\t%s\t%s\t%s\t%s\t0.00\t1.00\t1.00\tno\tno change\n" % (gene, count_dict[gene][0], count_dict[gene][1], fpkm_dict[gene][0], fpkm_dict[gene][1]))
        else:
            for gene in fpkm_dict.keys():
                w_text = "%s\t" % gene
                count_list = count_dict[gene][0] + count_dict[gene][1]
                fpkm_list = fpkm_dict[gene][0] + fpkm_dict[gene][1]
                for i in range(len(count_list)):
                    w_text += "%s\t%s\t" % (count_list[i], fpkm_list[i])
                w.write("%s%s\t%s\t0.00\t1.00\t1.00\tno\tno change\n" % (w_text, fpkm_mean[gene][0], fpkm_mean[gene][1]))

# stat_edger( 'genes.counts.matrix.P7_1_vs_P7_2.edgeR.DE_results','genes.counts.matrix', 'genes.TMM.EXPR.matrix', 'P1', 'P7', 'group.list.tmp', 0.05)
# stat_edger('/mnt/ilustre/users/bingxu.liu/project/2014peixun/denovo8/diffexpression/replicates_edgeR/unigene.counts.matrix.P1_vs_P7.edgeR.DE_results', '/mnt/ilustre/users/bingxu.liu/project/2014peixun/denovo8/diffexpression/unigene.counts.matrix', '/mnt/ilustre/users/bingxu.liu/project/2014peixun/denovo8/diffexpression/unigene.TMM.fpkm.matrix', 'P1', 'P7', '../group.list.tmp')
