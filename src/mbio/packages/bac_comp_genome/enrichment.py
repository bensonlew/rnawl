# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from scipy.stats import hypergeom
import pandas as pd
import statsmodels.stats.multitest as mt
import argparse
from table_format import annot_format
# from comp_function import *

methods = {'holm': 'h', 'hochberg': 'sh', 'hommel': 'ho',
           'bonferroni': 'b', 'bh': 'fdr_bh', 'by': 'fdr_by',
           'fdr': 'fdr_bh', 'none': 'none'}


def enrichment(tab4s, corrceted='fdr'):
    header = [
        'function', 'hitsInSet', 'numInSet', 'hitsInPop', 'numInPop', 'Pvalue'
    ]
    retval = []
    for i in tab4s.index:
        l = list(tab4s.loc[i])

        # hypergeom.sf(k, M, n, N)
        # k => hitsInSet - 1, M => numInPop, n => hitsInPop, N => numInSet
        if l[0] < 1:
            continue
        p = hypergeom.sf(l[0] - 1, l[1], l[2], l[3])
        retval.append([i, l[0], l[3], l[2], l[1], p])

    retval = pd.DataFrame(retval, columns=header).set_index('function')

    retval['Corrected'] = mt.multipletests(retval['Pvalue'],
                                           method=methods[corrceted])[1]

    return retval.sort_values(by='Corrected')


def format_tab(annot, tset, level):
    pop = []
    [pop.extend(ids.split(',')) for ids in annot['genelist']]
    pop = set(pop)  # 所有有功能注释的集合
    # print tset
    tset = tset & pop  # 用于富集的基因功能注释的集合
    # exit()
    pop_hits = {}  # key为功能类别, value为集合
    desc = {}
    desc_col = level + '_description'
    if level == 'Pathway':
        annot = annot.rename(columns={'Level3': desc_col})
    if desc_col not in annot.columns:
        annot[desc_col] = '-'
    for index in annot.index:
        # 注释文件中'genelist'列为每个功能下的基因列表，逗号隔开
        # desc_col 为待富集功能的功能描述
        i = annot.at[index, level]
        if i not in pop_hits:
            pop_hits[i] = set()
        pop_hits[i] |= set(annot.at[index, 'genelist'].split(','))
        desc[i] = annot.at[index, desc_col]
    tab4s = []
    col = [level, desc_col, 'genelist']
    enrich_info = []
    for k in pop_hits:
        tset_hits = set()
        tset_hits = tset & pop_hits[k]
        li = ','.join(tset_hits)
        tab4s.append([
            k, len(tset_hits), len(pop), len(pop_hits[k]), len(tset)
        ])
        enrich_info.append([
            k, desc[k], li
        ])

    tab4s = pd.DataFrame(data=tab4s).set_index(0)
    enrich_info = pd.DataFrame(data=enrich_info, columns=col).set_index(col[0])

    return tab4s, enrich_info


def one_background(annot, testsets, cor, level):
    out = None
    for tlabel, tset in testsets.items():
        tab4s, enrich_info = format_tab(annot, tset, level)
        enrich_result = enrichment(tab4s, cor)

        tmp_result = pd.concat([enrich_info, enrich_result], axis=1,
                               join='inner')
        if tlabel:
            tmp_result['label'] = tlabel
        out = pd.concat([out, tmp_result])
    out.index.name = level
    return out


def _main(args):
    testsets = {}  # 字典 key 为基因集名，value 基因列表
    with open(args.testset, 'r') as s:
        for l in s.readlines():
            line = l.strip().split()
            line.append(None)
            if line[1] not in testsets:
                testsets[line[1]] = set()
            testsets[line[1]].add(line[0])

    if args.formated == 'T':
        annot = pd.read_csv(args.annotable, sep='\t')
    else:
        annot = annot_format(args.annotable, args.functype)
        if args.sp != 'ALL':
            sp = args.sp.split(',')
            annot = annot[annot['sample'].isin(sp)]
    annot = annot[annot[args.level] != '-']

    bg_out = None
    if args.go != 'None':  # 类似GO注释有多个注释类别，分别作为背景进行富集
        grouped = annot.groupby(args.go)
        for tp, bg in grouped:
            tmp_out = one_background(bg, testsets, args.corrected, args.level)
            tmp_out[args.go] = bp
            bg_out = pd.concat([bg_out, tmp_out])
    else:
        bg_out = one_background(annot, testsets, args.corrected.lower(), args.level)

    bg_out.reset_index().to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-a', '--annotable', type=str, required=True,
                       help='功能注释table')
    parse.add_argument('-sp', default='ALL',
                       help='逗号分隔的字符串，对annotable的sample列筛选')
    parse.add_argument('-l', '--level', type=str, default=None, required=True,
                       help='选择要注释的功能层级所在列的列名')
    parse.add_argument('-s', '--testset', type=str, required=True,
                       help='用于富集分析的集合, 两列，第一列 基因名；\
                       第二列 集合名。如有表头，表头需以#开头')
    parse.add_argument('-o', '--out', type=str, default=None,
                       help='富集结果文件')
    parse.add_argument('-c', '--corrected', type=str, required=True,
                       help='矫正方法')
    parse.add_argument('-g', '--go', type=str, default='None',
                       help='功能注释表是否是类似GO具有多个注释类别的,\
                            若有则指出列名')
    parse.add_argument('-t', '--functype', type=str, default='KEGG',
                       help='功能的注释类长, KEGG, GO等')
    parse.add_argument('-f', '--formated', type=str, default='N',
                       help='注释表是否格式化，"N"表示未格式化，即第一列\
                            为基因名，注释信息在其他列；"Y"表示已经格式化，\
                            即第一列是注释信息，"genelist"列为逗号分隔的基因')

    args = parse.parse_args()
    _main(args)
