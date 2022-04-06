## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "HONGDONG"
# last_modify:20190305

import re
import argparse
import os
import gzip
from collections import defaultdict


def snp_type_distribution(infile, outfile_path):
    """
    sample_dict: {'AH03': {'CGGC': 286, 'CGAT': 445, 'TAAT': 488, 'TAGC': 405, 'CGTA': 1425, 'TACG': 1478},
    'AH19': {'CGGC': 274, 'CGAT': 397, 'TAAT': 470, 'TAGC': 389, 'CGTA': 1318, 'TACG': 1396},
    'CZ02': {'CGGC': 284, 'CGAT': 416, 'TAAT': 465, 'TAGC': 394, 'CGTA': 1420, 'TACG': 1445}}
    :param infile:
    :param outfile_path:
    :return:
    """
    with open(infile, 'r') as r:
        samples = []
        for line in r:
            if re.match('#CHROM.*', line):
                samples = line.strip().split('\t')
                break
        sample_dict = {}
        for i in range(9, len(samples)):
            # print samples[i]
            sample_dict[samples[i]] = {}
            sample_dict[samples[i]]['CGAT'] = 0
            sample_dict[samples[i]]['CGGC'] = 0
            sample_dict[samples[i]]['CGTA'] = 0
            sample_dict[samples[i]]['TAAT'] = 0
            sample_dict[samples[i]]['TACG'] = 0
            sample_dict[samples[i]]['TAGC'] = 0
        for line in r:
            if not re.match('#.*', line):
                temp = line.strip().split('\t')
                if temp[6] in ["PASS", "SNP", "INDEL", "FILTER"]:
                    type_ = ''.join([temp[3], temp[4]])
                    for i in range(9, len(temp)):
                        if temp[i].split(':')[0] not in ['0/0', './.', '*/*', '*/.']:
                            if type_ in ['CA', 'GT']:
                                sample_dict[samples[i]]['CGAT'] += 1
                            elif type_ in ['CG', 'GC']:
                                sample_dict[samples[i]]['CGGC'] += 1
                            elif type_ in ['CT', 'GA']:
                                sample_dict[samples[i]]['CGTA'] += 1
                            elif type_ in ['TA', 'AT']:
                                sample_dict[samples[i]]['TAAT'] += 1
                            elif type_ in ['TC', 'AG']:
                                sample_dict[samples[i]]['TACG'] += 1
                            elif type_ in ['TG', 'AC']:
                                sample_dict[samples[i]]['TAGC'] += 1
    outfile_name = os.path.join(outfile_path, 'snp_type_distribution.txt')
    with open(outfile_name, 'w') as w:
        w.write("#sample_name\tCG->AT\tCG->GC\tCG->TA\tTA->AT\tTA->CG\tTA->GC\n")
        for key in sample_dict.keys():
            w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(key, sample_dict[key]['CGAT'], sample_dict[key]['CGGC'],
                                                          sample_dict[key]['CGTA'], sample_dict[key]['TAAT'],
                                                          sample_dict[key]['TACG'], sample_dict[key]['TAGC']))


def indel_len_gene(infile, outfile_path):
    """
    sample_dict:
    {'AH03': {1: 76, 2: 35, 3: 24, 4: 6, 5: 6, 6: 5, 7: 3, 8: 3, 9: 1, 10: 0, -2: 26, -10: 1, -9: 1, -8: 3,
     -7: 0, -6: 6, -5: 6, -4: 14, -3: 16, -1: 91}, 'AH19': {1: 71, 2: 27, 3: 21, 4: 5, 5: 5, 6: 4, 7: 3, 8: 3, 9: 0,
     10: 0, -2: 23, -10: 2, -9: 1, -8: 3, -7: 0, -6: 6, -5: 5, -4: 12, -3: 17, -1: 82},
     'CZ02': {1: 72, 2: 35, 3: 23, 4: 7, 5: 6, 6: 4, 7: 2, 8: 3, 9: 1, 10: 0, -1: 84, -10: 1, -9: 1,
     -8: 4, -7: 0, -6: 4, -5: 4, -4: 9, -3: 19, -2: 26}}
    :param infile:
    :param outfile_path:
    :return:
    """
    with open(infile, 'r') as r:
        samples = []
        for line in r:
            if re.match('#CHROM.*', line):
                samples = line.strip().split('\t')
                break
        sample_dict = {}
        for i in range(9, len(samples)):
            # print samples[i]
            sample_dict[samples[i]] = {}
        for line in r:
            if not re.match('#.*', line):
                temp = line.strip().split('\t')
                if temp[6] in ["PASS", "SNP", "INDEL", "FILTER"]:
                    if re.search('intergenic_region|downstream_gene_variant|upstream_gene_variant', temp[7]):
                        continue
                    else:
                        indel_len = len(temp[3]) - len(temp[4])
                        for i in range(9, len(temp)):
                            if temp[i].split(':')[0] not in ['0/0', './.', '*/*', '*/.']:
                                if indel_len in sample_dict[samples[i]].keys():
                                    sample_dict[samples[i]][indel_len] += 1
                                else:
                                    sample_dict[samples[i]][indel_len] = 0
    outfile_name = os.path.join(outfile_path, 'indel_gene_distribution.txt')
    with open(outfile_name, 'w') as w:
        w.write("#sample_name\tlength\tindel_num\n")
        for key in sample_dict.keys():
            for key_ in sample_dict[key].keys():
                w.write("{}\t{}\t{}\n".format(key, key_, sample_dict[key][key_]))


if __name__ == '__main__':
    """
    对snp的突变类型进行统计，转换是同类碱基的置换（AT→GC及GC→AT）,颠换是不同类碱基的置换（AT→TA或CG,GC→CG或TA）
    CG->AT:['CA','GT']
    CG->GC:['CG','GC']
    CG->TA:['CT','GA']
    TA->AT:['TA','AT']
    TA->CG:['TC','AG']
    TA->GC:['TG','AC']
    """
    parser = argparse.ArgumentParser(description="用于统计snp的突变类型分布与indel基因区长度分布的统计")
    parser.add_argument("-i", "--infile", type=str, help="input vcf file:indel.anno.primary.vcf or snp.anno.primary.vcf")
    parser.add_argument("-t", "--type", type=str, help="snp or indel")
    parser.add_argument("-o", "--outfile", type=str, help="outfile path")
    args = parser.parse_args()
    if args.type == 'snp':
        snp_type_distribution(args.infile, args.outfile)
    elif args.type == 'indel':
        indel_len_gene(args.infile, args.outfile)
    else:
        raise Exception("-t 参数必须是snp或者indel")
