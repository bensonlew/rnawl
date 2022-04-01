#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import pandas as pd
import argparse
import json
import re
from table_format import *


def parse_corepan(annot, g2p, selection):
    tmp_colname = list(set(annot.columns) - set(['sample', 'genelist']))
    tmp_annot = []
    for index in annot.index:
        tmp_row = list(annot.loc[index, tmp_colname])
        sample = annot.at[index, 'sample']
        glist = annot.at[index, 'genelist']
        cat_set = {}
        for g in glist.split(','):
            cat = g2p[g]
            if cat not in cat_set:
                cat_set[cat] = set()
            cat_set[cat].add(g)
        for c in cat_set:
            if not selection:
                new_row = tmp_row + [sample + '-' + c, ','.join(cat_set[c])]
            elif c in selection:
                if len(selection) == 1:
                    new_row = tmp_row + [sample, ','.join(cat_set[c])]
                else:
                    new_row = tmp_row + [c, ','.join(cat_set[c])]
            else:
                continue
            tmp_annot.append(new_row)
    new_colname = tmp_colname + ['sample', 'genelist']
    annot = pd.DataFrame(data=tmp_annot, columns=new_colname)
    return annot


def gene2pancat(corepan, pancat):
    corepan = pd.read_csv(corepan, index_col=0, sep='\t')
    if set(['Sample_number', 'Gene_number']) < set(corepan.columns):
        corepan = corepan.drop(['Sample_number', 'Gene_number'], axis=1)
    if pancat == 'None':
        pancat = '{"uniqe": [1, 1], "dis": [1,99], "core": [100, 100]}'
    pancat = json.loads(pancat)
    max_v = max([float(v[1]) for v in pancat.values()])
    min_v = min([float(v[0]) for v in pancat.values()])
    print pancat
    allsamples = corepan.shape[1]
    g2p = {}
    for i in corepan.index:
        line = corepan.loc[i]
        s = line[line != '-']
        p = (len(s) + 0.0) / allsamples * 100
        if len(s) == 1:
            p = 1
        cat = which_cat(p, pancat, max_v, min_v)
        if not cat:
            exit('{} {}'.format(p, pancat))
        d = {g: cat for gl in list(s) for g in gl.split(',')}
        g2p.update(d)
    return g2p


def which_cat(p, pancat, max_v, min_v):
    for c, v in pancat.items():
        if max_v == float(v[1]) and p >= max_v:
            return c
        elif float(v[0]) <= p < float(v[1]):
            return c
        elif float(v[0]) == min_v and min_v >= p:
            return c


def corepan_sets(g2p):
    with open('corepan_sets.xls', 'w') as w:
        [w.write(g + '\t' + p + '\n') for g, p in g2p.iteritems()]


def get_temp_corepan(corepan):
    tmp = 'tmp_core_pan.xls'
    with open(corepan, 'r') as r, open(tmp, 'w') as w:
        repl = re.compile(r'[^\,\s]+\|')
        for l in r:
            new_l = repl.sub('', l)
            w.write(new_l)
    return tmp


def _main(args):
    if args.corepan:
        args.corepan = get_temp_corepan(args.corepan)
        g2p = gene2pancat(args.corepan, args.pancat)
        if args.corepansets:
            corepan_sets(g2p)
            return
        if args.selectedcat != 'ALL':
            selection = args.selectedcat.split(',')
        else:
            selection = None

    levels = args.level.lower().replace('cog', 'nog').split(',') if args.level else None
    annot = annot_format(args.annotable, args.functiontype, levels)
    print '###2###'
    print annot.columns
    sp_group = pd.read_csv(args.splist, header=None, sep='\t', comment='#')

    if args.groups.upper() != 'ALL':
        group_names = args.groups.split(',')
        sp_select = sp_group[sp_group[1].isin(group_names)]
        spl = list(sp_select[0])
    else:
        spl = list(sp_group[0])

    print '###3###'
    print set(annot['sample'])
    print spl
    annot = annot[annot['sample'].isin(spl)]
    print annot.columns
    print '###3###'

    if args.result_type in ['sum', 'mean', 'median']: # 样本名换成分组名
        annot = pd.merge(annot, sp_group, left_on='sample', right_on=0).drop(['sample', 0], axis=1).rename(columns={1: 'sample'})

    print '###4###'
    print annot.columns
    print len(set(annot['sample']))
    print '###4###'
    if args.corepan:
        annot = parse_corepan(annot, g2p, selection)
        if len(annot) == 0:
            raise Exception('当前选择的泛基因组分析结果基于选择的pangenome分类{}没有过滤得到genelist'.format(selection))

    annot = reshape(annot, 'sample', 'genelist', args.abundance, args.result_type)

    annot = annot.rename(columns={annot.columns[0]: 'Function'})
    annot.to_csv(args.output, index=None, sep='\t')


if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-a', '--annotable',
                       help='注释文件')
    parse.add_argument('-f', '--functiontype',
                       help='注释类别')
    parse.add_argument('-l', '--level', default=None,
                       help='选择要比较的基因注释层级，逗号分隔的字符串')
    parse.add_argument('-s', '--splist',
                       help='物种分组文件，第一列为物种，第二列为分组名')
    parse.add_argument('-sp', help='样本列表，逗号分隔的字符串')
    parse.add_argument('-g', '--groups', default='ALL', help='选择的分组名，逗号隔开')
    parse.add_argument('-c', '--corepan', default=None,
                       help='pangenome分析结果')
    parse.add_argument('-p', '--pancat', default='None',
                       help='pangenome分类方案')
    parse.add_argument('-o', '--output', default='out',
                       help='输出文件名')
    parse.add_argument('-b', '--abundance', default=None,
                       help='结果文件中输出以逗号分隔基因名的字符串\
                             还是基因的个数',)
    parse.add_argument('-d', '--selectedcat', default='ALL',
                       help='选择不同pangenome分类的结果')
    parse.add_argument('-r', '--result_type', default=None)
    parse.add_argument('--corepansets', default=None)
    args = parse.parse_args()

    _main(args)
