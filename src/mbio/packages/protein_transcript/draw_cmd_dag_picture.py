#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/19 12:30
@file    : draw_cmd_dag_picture.py
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

def draw_cpd_dag(source2target, node2color, out):
    # print(source2target)
    source2target_l = list()
    all_nodes = list()
    for s,t in source2target.items():
        all_nodes.append(s)
        if isinstance(t, list):
            for i in t:
                source2target_l.append((s,i))
                all_nodes.append(i)
        else:
            source2target_l.append((s,t))
            all_nodes.append(t)
    all_nodes = list(set(all_nodes))
    # print(all_nodes)
    G = pgv.AGraph(directed=True,strict=True)
    for s,t in source2target_l[::-1]:
        G.add_edge(#'\n'.join(textwrap.wrap(s, width=30, break_long_words=False)),
                   #'\n'.join(textwrap.wrap(t, width=30, break_long_words=False)),
                   t,
                   s,
                   label='', color='black',
                   style='solid',
                   minlen=1.5, arrowsize=1.3, penwidth=1.5)
    G.add_nodes_from(all_nodes)
    G.layout(prog='dot')
    G.graph_attr.update(dpi="100")
    G.node_attr.update(shape="box", style="rounded,filled", fillcolor="#FFFFFF")
    G.edge_attr.update(dir="back")
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

if __name__ == '__main__':

    if len(sys.argv) not in [4,5]:
        exit('USAGE:\npython %s relafile colorfile out_file' % sys.argv[0])
    relafile = sys.argv[1]
    colorfile = sys.argv[2]
    outfile = sys.argv[3]
    source2target = OrderedDict()
    node2color = dict()
    with open(relafile) as fr:
        _ = fr.readline()
        for line in fr:
            if line.strip():
                s,t = line.strip().split('\t')
                if not s in source2target:
                    source2target[s] = list()
                source2target[s].append(t)
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
    draw_cpd_dag(source2target, node2color, outfile)