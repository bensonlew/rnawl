# -*- coding: utf-8 -*-
import os, sys
import pandas as pd


def set_dict(dic, key, value):
    if key in dic:
        dic[key] += '###' + value
    else:
        dic[key] = value


def deal(cat, pro, out, d=True):
    desc = {}
    annot = {}
    out_list = []
    c_head = cat.readline().strip().split('\t')
    print c_head
    for line in cat:
        li = line.strip().split('\t')
        genes = li[2].split(',')
        map(lambda x: set_dict(annot, x, li[0]), genes)
        if d:
            #map(lambda x: set_dict(desc, x, li[-1]), genes)
            desc[li[0]] = li[-1]
    
    p_head = pro.readline().strip().split('\t')
    p_head[0] = c_head[0]
    print p_head
    for p in pro:
        pli = p.strip().split('\t')
        pli[1:] = map(int, pli[1:])
        if pli[0] in annot:
            pli[0] = annot[pli[0]]
            out_list.append(pli)

    print('pandas staff')
    tb = pd.DataFrame(data=out_list, columns=p_head)
    head0 = tb[p_head[0]].str.split('###', expand=True).stack().reset_index(level=1, drop=True).rename(p_head[0])
    tb = tb.drop(p_head[0], axis=1).join(head0)
    tb = tb.groupby(p_head[0]).agg(sum).reset_index()
    if d:
        tb['Description'] = tb[p_head[0]].apply(lambda x: desc[x])
    tb.to_csv(out, sep='\t', index=False)


def _main(cats, gpro, outs):
    for i in range(len(cats)):
        cat = open(cats[i], 'r')
        pro = open(gpro, 'r')
        if 'family_stat' in cats[i] or 'class_stat' in cats[i]:
            d = True
        else:
            d = False
        deal(cat, pro, outs[i], d)
        cat.close()
        pro.close()


if __name__ == '__main__':
    cats = sys.argv[1].split(',')
    gpro = sys.argv[2]
    outs = sys.argv[3].split(',')

    _main(cats, gpro, outs)
