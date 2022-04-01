#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/19 12:30
@file    : network_for_ppi.py
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
import argparse


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
    score_range = list(linspace(min(score_list),max(score_list), 1000))
    len_range = list(linspace(1, 3, 1000))
    score2len = zip(score_range, len_range)
    with open(rela_file) as rf:
        headers = rf.readline().strip().split('\t')
        infos = rf.readlines()
    for s,t in st_list:
        for line in infos:
            if not line.strip():
                continue
            if s+'\t' in line and t+'\t' in line:
                p, m, score = line.strip().split('\t')
                s = '\n'.join(textwrap.wrap(s, width=15, break_long_words=False))
                t = '\n'.join(textwrap.wrap(t, width=15, break_long_words=False))
                score = float(score)
                width = 1
                for n, c in enumerate(score2len):
                    if c[0] > score:
                        width = score2len[n-1][1]
                        break
                G.add_edge(
                    s,
                    t,
                    label='', color='black',
                    style='solid',
                    minlen=1.5, arrowsize=0, penwidth=width)
    G.add_nodes_from(node_list)
    G.layout(prog='dot')
    G.graph_attr.update(dpi="100")
    G.node_attr.update(shape="box", style="rounded,filled", fillcolor="#FFFFFF", height='3', width = '3')
    # G.node_attr.update(shape="circle", style="rounded,filled", fillcolor="#FFFFFF", height='3', width = '3')
    widths = list()
    for node in G.nodes():
        node = G.get_node(node)
        width = node.attr['width']
        widths.append(float(width))
    max_ = max(widths)
    edges = list(G.edges())
    for p in p2fc:
        if not p in G.nodes():
            fc_list.remove(p2fc[p])
    # diff_range = list(linspace(min(fc_list), 0.83, 500)) + list(linspace(1.2, max(fc_list), 500))
    down_range = list(linspace(min(fc_list), down, 500))
    up_range = list(linspace(up, max(fc_list), 500))
    # diff_range = list(linspace(min(fc_list), max(fc_list), 1000))
    # all_colors = list(Color.range_to(Color('green'), Color('red'), 1000))
    # all_colors = list(Color.range_to(Color('SpringGreen'), Color('DeepPink'), 1000))
    # all_colors = list(Color.range_to(Color('green'), Color('LightGreen'), 500)) + list(Color.range_to(Color('pink'), Color('red'), 500))
    down_colors = list(Color.range_to(Color('green'), Color('green', luminance=0.9), 500))
    up_colors = list(Color.range_to(Color('red', luminance=0.9), Color('red'), 500))
    upzip = zip(up_range,up_colors)
    downzip = zip(down_range,down_colors)
    # range2color = zip(diff_range, all_colors)
    for node in G.nodes():
        # print(node)
        count = 0
        for e in edges:
            if node in e:
                count += 1
        fc = 0
        if node in p2fc:
            fc = p2fc[node]
            # print(fc)
        color = 'white'
        if fc and fc <= down:
            for n, c in enumerate(downzip):
                if c[0] >= fc:
                    color = downzip[n - 1][1].hex_l
                    break
        if fc and fc >= up:
            for n, c in enumerate(upzip):
                if c[0] >= fc:
                    color = upzip[n - 1][1].hex_l
                    break
        # print(color)
        node = G.get_node(node)
        # node.attr['shape'] = 'circle'
        node.attr['shape'] = 'box'
        node.attr['fillcolor'] = color
        # node.attr['width'] = str(float(node.attr['width']) + math.log(count, 10000))
        node.attr['width'] = str(max_ + math.log(count, 10))
        node.attr['height'] = str(float(node.attr['height'])/3 + math.log(count, 10))
        # node.attr['height'] = node.attr['width']

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
            if (p,m) in pm_list or (m, p) in pm_list:
                continue
            pm_list.append((p,m))
            for nt, tree in enumerate(tree_list):
                # print(tree)
                if tree[0] == m:
                    flag = 1
                    tree_list_c[nt] = [p] + tree
                elif tree[0] == p:
                    flag = 1
                    tree_list_c[nt] = [m] + tree
                elif tree[-1] == p:
                    flag = 1
                    tree_list_c[nt] = tree + [m]
                elif tree[-1] == m:
                    flag = 1
                    tree_list_c[nt] = tree + [p]
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
    # print(tree_list)
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
                    st = (tree[n], tree[n + 1])
                    ts = (tree[n + 1], tree[n])
                    if st not in st_list and ts not in st_list:
                        st_list.append(st)
                    if s not in node_list:
                        node_list.append(s)
                except:
                    if tree[-1] not in node_list:
                        node_list.append(tree[-1])
    # max_ = max([len(tree) for tree in tree_list])
    # twos = list()
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

def predeal_ppi(ppi_file):
    ppi_pd = pd.read_csv(ppi_file, sep='\t')
    while ppi_pd.shape[0] > 250:
        min_ = ppi_pd.iloc[:,2].min()
        ppi_pd = ppi_pd[ppi_pd.iloc[:, 2] > (min_+1)]
    new_ppi = ppi_file + '_for_picture.xls'
    ppi_pd.to_csv(new_ppi,sep='\t',index=False,header=True)
    return new_ppi


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="generate ppi picture")
    parser.add_argument("-ppi_info", type=str, required=True, help="ppi info file")
    parser.add_argument("-de_info", type=str, required=True,
                        help="deinfo file")
    parser.add_argument("-up", type=float, default=1.2, help="the up value")
    parser.add_argument("-down", type=float, default=0.83, help='the down value')
    parser.add_argument("-out", type=str, default='ppi',help='the prefix of out file')
    parser.add_argument("-draw_two", type=str, default='no',
                        help='whether to draw two nodes')

    args = parser.parse_args()
    ppifile = args.ppi_info
    deinfo = args.de_info
    outfile = args.out
    up = args.up
    down = args.down
    draw_two = args.draw_two
    ppifile = predeal_ppi(ppifile)
    de_pd = pd.read_csv(deinfo, sep='\t')
    fc_list = de_pd.iloc[:,1].to_list()
    p2fc = {z[0]: z[1] for z in zip(de_pd.iloc[:, 0], de_pd.iloc[:, 1])}
    ppi_pd = pd.read_csv(ppifile, sep='\t')
    score_list = ppi_pd.iloc[:,2].to_list()
    st_list, node_list = analyse_corr(ppifile)
    # print(st_list, node_list)
    draw_network(ppifile, st_list, node_list, outfile)