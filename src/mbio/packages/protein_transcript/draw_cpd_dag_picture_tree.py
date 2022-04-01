#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/19 12:30
@file    : draw_cpd_dag_picture_tree.py
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
from collections import OrderedDict
from numpy import linspace


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

def draw_cpd_dag(st_list, node_list, node2color, out):

    G = pgv.AGraph(directed=True,strict=True)
    for s,t in st_list:
        G.add_edge(#'\n'.join(textwrap.wrap(s, width=30, break_long_words=False)),
                   #'\n'.join(textwrap.wrap(t, width=30, break_long_words=False)),
                   s,
                   t,
                   label='', color='black',
                   style='solid',
                   minlen=1.5, arrowsize=1.3, penwidth=1.5)
    G.add_nodes_from(node_list)
    G.layout(prog='dot')
    G.graph_attr.update(dpi="100")
    G.node_attr.update(shape="box", style="rounded,filled", fillcolor="#FFFFFF")
    # G.edge_attr.update(dir="back")
    for n,c in node2color.items():
        color = 'white'
        # if float(c[1]) < 0.9:
        #     v = float(c[1]) + 0.1
        # else:
        #     v = float(c[1])
        v = float(c[1])
        if c[0] == 'blue':
            for n_,tur in enumerate(blue2range):
                if tur[0] > v:
                    color = blue2range[n_-1][1].hex_l
                    break
        else:
            for n_,tur in enumerate(red2range):
                if tur[0] > v:
                    color = red2range[n_-1][1].hex_l
                    break
        # print(color)
        try:
            node = G.get_node(n)
            node.attr['fillcolor'] = color
        except:
            print(n+'不在node里面，为啥还有颜色')
    G.draw(out + '.tmp.svg', prog="dot")
    G.draw(out + '.pdf', prog="dot")
    G.draw(out + '.png', prog="dot")
    moidify_svg(out + '.tmp.svg', out + '.svg')
    os.system("rm " + out + '.tmp.svg')

def analyse_tree(rela_file):
    st_list = list()
    tree_list = list()
    node_list = list()
    with open(rela_file) as rf:
        _ = rf.readline()
        for line in rf:
            if line.strip():
                tree_list_c = tree_list[::]
                s, t = line.strip().split('\t')
                flag = 0
                for nt, tree in enumerate(tree_list):
                    for n, node in enumerate(tree):
                        if node == t:
                            flag = 1
                            if n != 0:
                                if tree[n - 1] == s:
                                    flag = 1
                                    continue
                                else:
                                    tree_list_c.append([s] + tree[n:])
                            else:
                                tree_list_c[nt] = [s] + tree
                    if tree[-1] == s:
                        tree_list_c[nt] += [t]
                if not flag:
                    tree_list_c.append([s, t])
                tree_list = tree_list_c[::]
    # print(tree_list)
    max_ = max([len(tree) for tree in tree_list])
    for n in range(max_):
        for tree in tree_list:
            try:
                s = tree[n]
                st = (tree[n], tree[n + 1])
                if st not in st_list:
                    st_list.append(st)
                if s not in node_list:
                    node_list.append(s)
            except:
                if tree[-1] not in node_list:
                    node_list.append(tree[-1])

    return st_list, node_list

if __name__ == '__main__':

    if len(sys.argv) not in [4,5]:
        exit('USAGE:\npython %s relafile colorfile out_file' % sys.argv[0])
    relafile = sys.argv[1]
    colorfile = sys.argv[2]
    outfile = sys.argv[3]
    source2target = OrderedDict()
    node2color = dict()
    st_list, node_list = analyse_tree(relafile)
    with open(colorfile) as cr:
        _ = cr.readline()
        for line in cr:
            if line.strip():
                n,c,v = line.strip().split('\t')
                node2color[n] = [c,v]
    _range = list(linspace(0, 1, 1000))
    all_red_colors = list(Color.range_to(Color('Pink'), Color('red'), 1000))
    all_blue_colors = list(Color.range_to(Color('PowderBlue'), Color('blue'), 1000))
    red2range = zip(_range, all_red_colors)
    blue2range = zip(_range, all_blue_colors)
    draw_cpd_dag(st_list, node_list, node2color, outfile)