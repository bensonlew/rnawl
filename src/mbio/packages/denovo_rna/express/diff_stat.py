# !/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "qiuping"
# last_modify:20161101
from collections import defaultdict
import math


class Express(object):
    def __init__(self):
        self.gene = None
        self.counts = {}  # 每个样本对应的count值
        self.fpkms = {}  # 每个样本对应的fpkm值


class DiffStat(object):
    def __init__(self):
        self.express_info = {}
        self.samples = []
        self.group_info = defaultdict(list)

    def get_express_info(self, countfile, fpkmfile):
        """
        countfile：基因计数矩阵表
        fpkmfile：基因表达量矩阵表
        return：express_info字典，键为基因名称，值为Express对象
        """
        with open(countfile, 'rb') as c, open(fpkmfile, 'rb') as f:
            head = c.readline().strip('\n').split('\t')
            self.samples = head[1:]
            for line in c:
                line = line.strip('\n').split('\t')
                gene = line[0]
                foo = Express()
                foo.gene = gene
                for sam in self.samples:
                    foo.counts[sam] = float(line[head.index(sam)])
                self.express_info[gene] = foo

            f.readline()
            for line in f:
                line = line.strip('\n').split('\t')
                gene = line[0]
                tmp = {}
                for sam in self.samples:
                    tmp[sam] = float(line[head.index(sam)])
                self.express_info[gene].fpkms = tmp
        return self.express_info

    def get_group_info(self, group):
        """
        group:分组文件
        return：字典，键为分组名称，值为该分组对应的所有样本名的列表；
        eg: {'A':[1,2,3], 'B':[4,5,6]}
        """
        with open(group, 'rb') as g:
            g.readline()
            for line in g:
                line = line.strip('\n').split()
                self.group_info[line[1]].append(line[0])
        return self.group_info

    def get_mean(self, samples, count_dict, fpkm_dict):
        """
        samples: 样本的列表
        count_dict：全部样本对应的counts值的字典
        fpkm_dict：全部样本对应的fpkm值的字典
        return：samples中的样本的count、fpkm的均值
        """
        sum_count = 0
        sum_fpkm = 0
        for sam in samples:
            sum_count += count_dict[sam]
            sum_fpkm += fpkm_dict[sam]
        mean_count = float(sum_count) / len(samples)
        mean_fpkm = float(sum_fpkm) / len(samples)
        return round(mean_count, 3), round(mean_fpkm, 3)

    def diff_stat(self, express_info, edgr_result, control, other, output, group_info=None, regulate=True, diff_ci=0.05, fc=0):
        """
        express_info:字典，键为基因名，值为Express对象
        edgr_result:edgr分析得到的结果文件
        control:对照组/样本名
        other:实验组/样本名
        group_info:字典，键为分组名称，值为该分组对应的所有样本名的列表
        regulate:是否做上下调分析
        diff_ci:显著差异水平
        fc:差异倍数，当logfc大于该值时则认为上调，否则下调
        """
        with open(edgr_result, 'rb') as r, open('%s/%s_vs_%s_edgr_stat.xls' % (output, control, other), 'wb') as w:
            r.readline()
            if group_info:
                head = "gene_id\t%s_mean_count\t%s_mean_count\t%s_mean_fpkm\t%s_mean_fpkm\tlog2fc(%s/%s)\tlogCPM(%s*%s)\tpvalue\tfdr\tsignificant\tregulate\n" % (control, other, control, other, other, control,other, control)
            else:
                head = "gene_id\t%s_count\t%s_count\t%s_fpkm\t%s_fpkm\tlog2fc(%s/%s)\tlogCPM(%s*%s)\tpvalue\tfdr\tsignificant\tregulate\n" % (control, other, control, other, other, control,other, control)
            w.write(head)
            for line in r:
                line = line.strip('\n').split('\t')
                gene = line[0]
                pvalue = float(line[-2])
                logCPM = float(line[-3])
                fdr = float(line[-1])
                counts = express_info[gene].counts
                fpkms = express_info[gene].fpkms
                if group_info:
                    con_sams = group_info[control]
                    oth_sams = group_info[other]
                    control_count, control_fpkm = self.get_mean(con_sams, counts, fpkms)
                    other_count, other_fpkm = self.get_mean(oth_sams, counts, fpkms)
                else:
                    control_count = express_info[gene].counts[control]
                    control_fpkm = express_info[gene].fpkms[control]
                    other_count = express_info[gene].counts[other]
                    other_fpkm = express_info[gene].fpkms[other]
                #lfc = (other_fpkm + 0.1) / (control_fpkm + 0.1)
                #logfc = round(math.log(lfc, 2), 3)
                logfc = float(line[-4])
                if regulate:
                    if logfc > fc:
                        reg = 'up'
                    elif logfc < fc:
                        reg = 'down'
                    else:
                        reg = 'no change'
                else:
                    reg = 'undone'

                if fdr < diff_ci:
                    sig = 'yes'
                else:
                    sig = 'no'
                w.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene, control_count, other_count, control_fpkm, other_fpkm, '%0.4g' %logfc,'%0.4g' %logCMP, '%0.4g' % pvalue, '%0.4g' % fdr, sig, reg))

# a = DiffStat()
# a.get_express_info('/mnt/ilustre/users/sanger-dev/workspace/20161101/TestBase_tsn_50/ExpAnalysis/MergeRsem/genes.counts.matrix', '/mnt/ilustre/users/sanger-dev/workspace/20161101/TestBase_tsn_50/ExpAnalysis/MergeRsem/genes.TMM.fpkm.matrix')
# # print a.express_info
# a.diff_stat(express_info=a.express_info, edgr_result='/mnt/ilustre/users/sanger-dev/workspace/20161101/TestBase_tsn_50/ExpAnalysis/DiffExp/edger_result/genes.counts.matrix.E18_1_vs_E18_2.edgeR.DE_results', control='E18_1', other='E18_2', output='/mnt/ilustre/users/sanger-dev/workspace/20161101/TestBase_tsn_50/ExpAnalysis/MergeRsem/', group_info=None, regulate=True, diff_ci=0.05)
