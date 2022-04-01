#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/19 12:30
@file    : network_for_express.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import sys
sys.path.append('/mnt/ilustre/centos7users/yitong.feng/.pythonlib/lib/python2.7/site-packages/')
import pygraphviz as pgv
from colour import Color
import textwrap
import os
import xml.etree.ElementTree as ET
from numpy import linspace
import pandas as pd
from sklearn import preprocessing
import math


def moidify_svg(svg_old, svg_new):
    """
    修改svg展示，防止图片超出边界 liubinxu
    """
    xml = ET.parse(svg_old)
    root = xml.getroot()
    root.attrib['viewBox'] = '0 0 {} {}'.format(root.attrib['width'][:-2], root.attrib['height'][:-1])
    g1 = root.find("{http://www.w3.org/2000/svg}g")
    gsubs = g1.findall("{http://www.w3.org/2000/svg}g")
    for gsub in gsubs:
        texts= gsub.findall('{http://www.w3.org/2000/svg}text')
        for text in texts:
            text.attrib['font-size'] = str(float(text.attrib['font-size']) - 2)
    xml.write(svg_new)

def draw_network(rela_file, st_list, node_list, out):
    node_list = ['\n'.join(textwrap.wrap(n, width=15, break_long_words=False)) for n in node_list]
    G = pgv.AGraph()
    p_list = list()
    m_list = list()
    _range = list(linspace(-1, 1, 1000))
    all_colors = list(Color.range_to(Color('SpringGreen'), Color('DeepPink'), 1000))
    range2color = zip(_range, all_colors)
    with open(rela_file) as rf:
        headers = rf.readline().strip().split('\t')
        infos = rf.readlines()
    for s,t in st_list:
        for line in infos:
            if not line.strip():
                continue
            if s+'\t' in line and t+'\t' in line:
                p, m, cor, pv, fdr = line.strip().split('\t')
                p = '\n'.join(textwrap.wrap(p, width=15, break_long_words=False))
                m = '\n'.join(textwrap.wrap(m, width=15, break_long_words=False))
                s = '\n'.join(textwrap.wrap(s, width=15, break_long_words=False))
                t = '\n'.join(textwrap.wrap(t, width=15, break_long_words=False))
                p_list.append(p)
                m_list.append(m)
                cor = float(cor)
                color = 'black'
                for n, c in enumerate(range2color):
                    if c[0] > cor:
                        color = range2color[n-1][1].hex_l
                        break
                G.add_edge(
                    s,
                    t,
                    label='', color=color,
                    style='solid',
                    minlen=1.5, arrowsize=0, penwidth=2.5)
    G.add_nodes_from(node_list)
    G.layout(prog='dot')
    G.graph_attr.update(dpi="100")
    G.node_attr.update(shape="box", style="rounded,filled", fillcolor="#FFFFFF", height='3', width = '3')
    widths = list()
    for node in G.nodes():
        node = G.get_node(node)
        width = node.attr['width']
        widths.append(float(width))
    max_ = max(widths)
    edges = list(G.edges())
    for p in p_list:
        count = 0
        for e in edges:
            if p in e:
                count += 1
        node = G.get_node(p)
        node.attr['shape'] = 'triangle'
        node.attr['fillcolor'] = Color('Cyan').hex_l
        # node.attr['width'] = str(float(node.attr['width']) + math.log(count, 10000))
        node.attr['width'] = str(max_ + math.log(count, 10))
        # node.attr['height'] = str(float(node.attr['height']) + math.log(count, 10))
        node.attr['height'] = node.attr['width']
    for m in m_list:
        count = 0
        for e in edges:
            if m in e:
                count += 1
        node = G.get_node(m)
        node.attr['shape'] = 'circle'
        node.attr['fillcolor'] = Color('Pink').hex_l
        # node.attr['width'] = str(float(node.attr['width'])+math.log(count, 10000))
        node.attr['width'] = str(max_ + math.log(count, 10))
        # node.attr['height'] = str(float(node.attr['height'])+math.log(count, 10))
        node.attr['height'] = node.attr['width']
    G.draw(out + '.tmp.svg', prog="dot")
    G.draw(out + '.pdf', prog="dot")
    G.draw(out + '.png', prog="dot")
    moidify_svg(out + '.tmp.svg', out + '.svg')
    os.system("rm " + out + '.tmp.svg')

def analyse_corr(rela_file):
    st_list = list()
    tree_list = list()
    node_list = list()
    pm_list = list()
    with open(rela_file) as rf:
        _ = rf.readline()
        for line in rf:
            if not line.strip():
                continue
            flag = 0
            tree_list_c = tree_list[::]
            p, m = line.strip().split('\t')[:2]
            # print(p,m)
            if (p,m) in pm_list:
                continue
            pm_list.append((p,m))
            for nt, tree in enumerate(tree_list):
                # print(tree)
                if tree[0] == m:
                    flag = 1
                    tree_list_c[nt] = [p] + tree
                    # continue
                elif tree[0] == p:
                    flag = 1
                    tree_list_c[nt] = [m] + tree
                    # continue
                elif tree[-1] == p:
                    flag = 1
                    tree_list_c[nt] = tree + [m]
                    # continue
                elif tree[-1] == m:
                    flag = 1
                    tree_list_c[nt] = tree + [p]
                    # continue
                else:
                    for n, node in enumerate(tree):
                        if node == p:
                            tmp = tree[:n+1] + [m]
                            if tmp not in tree_list_c:
                                tree_list_c.append(tmp)
                            flag = 1
                            break
                        if node == m:
                            tmp = tree[:n+1] + [p]
                            if tmp not in tree_list_c:
                                tree_list_c.append(tmp)
                            flag = 1
                            break
            if not flag:
                tree_list_c.append([p,m])
            tree_list = tree_list_c[::]
            # tree_list = list(set(tree_list_c))
    # print(tree_list)
    # max_ = max([len(tree) for tree in tree_list])
    len2tree = dict()
    for tree in tree_list:
        if len(tree) not in len2tree:
            len2tree[len(tree)] = list()
        len2tree[len(tree)].append(tree)
    twos = list()
    if 2 in len2tree:
        twos = len2tree[2]
        del len2tree[2]
    for l in sorted(len2tree.keys(), reverse=True):
        for tree in len2tree[l]:
            for n in range(l):
                try:
                    s = tree[n]
                    st = (tree[n], tree[n+1])
                    ts = (tree[n+1], tree[n])
                    if st not in st_list and ts not in st_list:
                        st_list.append(st)
                    if s not in node_list:
                        node_list.append(s)
                except:
                    if tree[-1] not in node_list:
                        node_list.append(tree[-1])
    # for n in range(max_):
    #     for tree in tree_list:
    #         if len(tree) == 2:
    #             if tree not in twos:
    #                 twos.append(tree)
    #             continue
    #         try:
    #             s = tree[n]
    #             st = (tree[n], tree[n+1])
    #             ts = (tree[n+1], tree[n])
    #             if st not in st_list and ts not in st_list:
    #                 st_list.append(st)
    #             if s not in node_list:
    #                 node_list.append(s)
    #         except:
    #             if tree[-1] not in node_list:
    #                 node_list.append(tree[-1])
    if draw_two.lower() == 'yes':
        for s,t in twos:
            if (s,t) not in st_list and (t,s) not in st_list:
                st_list.append((s,t))
            if s not in node_list:
                node_list.append(s)
            if t not in node_list:
                node_list.append(t)

    return st_list, node_list

def predeal_relafile(rela_file):
    rela_pd = pd.read_csv(rela_file, sep='\t')
    while rela_pd.shape[0] > 250:
        pos_min = min([x for x in rela_pd.iloc[:,2] if x > 0])
        neg_max = max([x for x in rela_pd.iloc[:,2] if x < 0])
        rela_pd = rela_pd[((rela_pd.iloc[:, 2] > (pos_min+0.001)) | (rela_pd.iloc[:, 2] < (neg_max-0.001)))]
    new_rela = rela_file + '_for_picture.xls'
    rela_pd.to_csv(new_rela,sep='\t',index=False,header=True)
    return new_rela


if __name__ == '__main__':

    if len(sys.argv) not in [3, 4]:
        exit('USAGE:\npython %s relafile out_file' % sys.argv[0])
    relafile = sys.argv[1]
    relafile = predeal_relafile(relafile)
    outfile = sys.argv[2]
    try:
        draw_two = sys.argv[3]
    except:
        draw_two = 'yes'
    st_list, node_list = analyse_corr(relafile)
    # print(st_list, node_list)
    draw_network(relafile, st_list, node_list, outfile)