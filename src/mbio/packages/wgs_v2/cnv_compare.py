## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "HONGDONG"
# last_modify:20190307

import re
import argparse
import os
from collections import defaultdict


class CnvCompare(object):
    def __init__(self, sample1_cnv, sample2_cnv, region, marktype, output_dir, region_type='real_region'):
        """
        :param sample1_cnv: 样本1的cnv.anno.xls
        :param sample2_cnv: 样本2的cnv.anno.xls
        :param region: 基因组区域
        :param marktype: 样本间是去相同还是不同还是都选
        :param output_dir: 输出结果目录
        :param region_type: 区域筛选方式，real_region时chr1:1-2是指一个范围，如果不是则指具体一个snp或者indel
        """
        self.check_file_exitst(sample1_cnv)
        self.check_file_exitst(sample2_cnv)
        self.sample1_cnv = sample1_cnv
        self.sample2_cnv = sample2_cnv
        self.region = region
        self.marktype = marktype
        self.output_dir = output_dir
        self.region_type = region_type
        self.name1 = os.path.basename(self.sample1_cnv).split('.')[0]
        self.name2 = os.path.basename(self.sample2_cnv).split('.')[0]

    def change_region(self, chr, start, end, type_='real_region'):
        """
        "chr1:1-2,chr2:2-4"
        :return: True 该位点在范围内
        """
        if self.region.lower() == 'all':
            return True
        else:
            for reg in self.region.split(','):
                if type_ == 'real_region':
                    try:
                        start_ = int(reg.split(':')[1].split('-')[0])
                    except:
                        start_ = 0
                    try:
                        end_ = int(reg.split(':')[1].split('-')[1])
                    except:
                        end_ = 100000000000000000000000000000000000000000000
                    if chr == reg.split(':')[0] and int(start) >= start_ and int(end) <= end_:
                        return True
                else:
                    if chr + ':' + str(start) + '-' + str(end) == reg:
                        return True
        return False

    def read_sample(self, sample_cnv, name):
        cnv = defaultdict(list)
        with open(sample_cnv, 'r') as r:
            for line in r:
                if not re.match('#.*', line):
                    temp = line.strip().split('\t')
                    gene_detail = '--' if temp[6] == '0' else temp[7]
                    if self.change_region(temp[0], temp[1], temp[2], self.region_type):
                        cnv['\t'.join([temp[0], temp[1], temp[2], temp[3], temp[4]])] = \
                            [name, temp[5], temp[6], gene_detail]
        return cnv

    def cnv_compare(self):
        outfile = os.path.join(self.output_dir, '{}_VS_{}_{}.xls'.format(self.name1, self.name2, self.marktype))
        with open(outfile, 'w') as w:
            w.write('#chr\tstart\tend\tlength\ttype\t{}_genotype\t{}_genotype\t{}_pvalue\t{}_pvalue\tgene_num\tgene\n'
                    .format(self.name1, self.name2, self.name1, self.name2))
            cnv1 = self.read_sample(self.sample1_cnv, self.name1)
            cnv2 = self.read_sample(self.sample2_cnv, self.name2)
            if self.marktype == 'same':  # 样本1与样本2相同的位点
                for key1 in sorted(cnv1.keys()):
                    if key1 in cnv2.keys():
                        w.write("{}\tT\tT\t{}\t{}\t{}\t{}\n".format(key1, cnv1[key1][1], cnv2[key1][1], cnv1[key1][2],
                                                                    cnv1[key1][3]))
            elif self.marktype == 'diff':
                for key1 in sorted(cnv1.keys()):   # 样本1在样本2中不同的点
                    if key1 not in cnv2.keys():
                        w.write("{}\tT\tF\t{}\t--\t{}\t{}\n".format(key1, cnv1[key1][1], cnv1[key1][2], cnv1[key1][3]))
                for key1 in sorted(cnv2.keys()):   # 样本2在样本1中不同的点
                    if key1 not in cnv1.keys():
                        w.write("{}\tF\tT\t--\t{}\t{}\t{}\n".format(key1, cnv2[key1][1], cnv2[key1][2], cnv2[key1][3]))
            else:
                for key1 in sorted(cnv1.keys()):
                    if key1 in cnv2.keys():
                        w.write("{}\tT\tT\t{}\t{}\t{}\t{}\n".format(key1, cnv1[key1][1], cnv2[key1][1], cnv1[key1][2],
                                                                    cnv1[key1][3]))
                for key1 in sorted(cnv1.keys()):
                    if key1 not in cnv2.keys():
                        w.write("{}\tT\tF\t{}\t--\t{}\n".format(key1, cnv1[key1][1], cnv1[key1][2], cnv1[key1][3]))
                for key1 in sorted(cnv2.keys()):
                    if key1 not in cnv1.keys():
                        w.write("{}\tF\tT\t--\t{}\t{}\n".format(key1, cnv2[key1][1], cnv2[key1][2], cnv2[key1][3]))
        self.cnv_stat(outfile)

    def check_file_exitst(self, file_path):
        if not os.path.exists(file_path):
            raise Exception("文件{}不存在！".format(file_path))

    def cnv_stat(self, file_path):
        """
        对cnv结果进行差异统计
        :return:
        """
        chr_info = defaultdict(list)
        out = file_path.split('.')[0] + '.stat.xls'
        with open(file_path, 'r') as r, open(out, 'w') as w:
            w.write("#chrid\tdeletion\tduplication\tgene\n")
            for line in r:
                if not re.match('^#.*', line):
                    temp = line.strip().split('\t')
                    if temp[0] in chr_info.keys():
                        if temp[4] == 'deletion':
                            chr_info[temp[0]][0] += 1
                        if temp[4] == 'duplication':
                            chr_info[temp[0]][1] += 1
                        # noinspection PyBroadException
                        try:
                            chr_info[temp[0]][2] += int(temp[9])
                        except:
                            chr_info[temp[0]][2] += 0
                    else:
                        # noinspection PyBroadException
                        try:
                            chr_info[temp[0]] = [1, 1, int(temp[9])]
                        except:
                            chr_info[temp[0]] = [1, 1, 0]
            # for key in sorted(chr_info.keys()):
            for key in self.sort_chr_sca(chr_info):
                w.write('{}\t{}\t{}\t{}\n'.format(key, chr_info[key][0], chr_info[key][1], chr_info[key][2]))

    def sort_chr_sca(self, chr_info):
        chrs = []
        scas = []
        other = []
        for key in chr_info.keys():
            if key.startswith('chr'):
                chrs.append(key)
            elif key.startswith('sca'):
                scas.append(key)
            else:
                other.append(key)
        chrs = self.sort_list(chrs)
        scas = self.sort_list(scas)
        if other:
            other = self.sort_list(other)
        chrs.extend(scas)
        chrs.extend(other)
        return chrs

    def sort_list(self, list_):
        # noinspection PyBroadException
        try:
            list_.sort(key=lambda i: int(re.findall("\d+", i)[0]))
        except:
            list_.sort()
        return list_


if __name__ == '__main__':
    """
    用于两个样本的cnv比较分析统计
    python cnv_compare.py -i ./ -r all -m 'same,diff' -s 'AH03|AH19,CZ02|AH03' -o ./
    """
    parser = argparse.ArgumentParser(description="用于cnv进行样本间的比较分析example python cnv_compare.py -i ./ -r all"
                                                 " -m 'same,diff' -s 'AH03|AH19,CZ02|AH03' -o ./")
    parser.add_argument("-i", "--infile", type=str, help="input cnv anno xls file path！")
    parser.add_argument("-r", "--region", type=str, help="select genome region, chr1:1-500,chr2:2-4")
    parser.add_argument("-m", "--marktype", type=str, help="same or diff or all, same,diff")
    parser.add_argument("-s", "--samples", type=str, help="a|b,b|c")
    parser.add_argument("-rt", "--region_type", type=str, help="real_region or not real region", default='real_region')
    parser.add_argument("-o", "--outfile", type=str, help="outfile path")
    args = parser.parse_args()
    samples = args.samples.split(',')
    marktype = args.marktype.split(',')
    if len(samples) != len(marktype):
        raise Exception("样本对数与比较类型个数不一致")
    for i in range(0, len(samples)):
        sample1, sample2 = samples[i].split('|')
        sample1_path = os.path.join(args.infile, '{}.cnv.anno.xls'.format(sample1))
        sample2_path = os.path.join(args.infile, '{}.cnv.anno.xls'.format(sample2))
        a = CnvCompare(sample1_path, sample2_path, args.region, marktype[i], args.outfile, args.region_type)
        a.cnv_compare()
