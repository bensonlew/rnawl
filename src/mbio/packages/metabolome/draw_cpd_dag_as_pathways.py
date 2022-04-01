# -*- coding: utf-8 -*-
"""
@time    : 2019/4/19 12:30
@file    : draw_cpd_dag_picture_tree.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import sys
#sys.path.append('/mnt/ilustre/centos7users/yitong.feng/.pythonlib/lib/python2.7/site-packages/')
import pygraphviz as pgv
from colour import Color
import textwrap
import os
import xml.etree.ElementTree as ET
from collections import OrderedDict
from numpy import linspace
import json
import re
import shutil


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
    for s, t in st_list:
        G.add_edge(#'\n'.join(textwrap.wrap(s, width=30, break_long_words=False)),
                   #'\n'.join(textwrap.wrap(t, width=30, break_long_words=False)),
                   s,
                   t,
                   label='', color='black',
                   style='solid',
                   minlen=1.5, arrowsize=1.3, penwidth=1.5)

    G.add_nodes_from(node_list)
    ###G.layout(prog='dot')
    G.graph_attr.update(dpi="100")
    G.node_attr.update(shape="box", style="rounded,filled", fillcolor="#FFFFFF")
    for n,c in node2color.items():
        color = 'red'
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
    return


def analyse_tree(net):
    st_list = list()
    tree_list = list()
    node_list = list()
    analysed = list()
    for rela in net:
        tree_list_c = tree_list[::]
        s, t = rela[0], rela[1]
        s = re.sub(r'cpd:', '', s)
        t = re.sub(r'cpd:', '', t)
        flag = 0
        if (s, t) in analysed:
            continue
        else:
            analysed.append((s, t))
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
    os.environ['LD_LIBRARY_PATH'] = "/usr/lib/:/usr/lib64:" + os.environ['LD_LIBRARY_PATH']
    cwd = os.getcwd()
    #db_dir = '/mnt/ilustre/users/ting.kuang/ITRAQ/db/KEGG/'
    #db_dir = '/mnt/ilustre/users/sanger-dev/app/database/metabolome/topo_json/'
    db_dir = sys.argv[3]    ## database path
    #org = '2019-01.ko'
    #org = 'hsa'
    org = sys.argv[4]  # org name , is subfolder of db_dir
    colorfile = sys.argv[1]
    outpath = sys.argv[2]
    if len(sys.argv) == 6:
        org_file = sys.argv[5]
    else:
        org_file = ''

    if os.path.exists(outpath):
        shutil.rmtree(outpath)
    os.mkdir(outpath)
    primary_list = ['01100', '01110', '01120', '01130', '01200', '01210', '01212', '01230', '01220']
    with open(db_dir + org + "/network.json", 'r') as load_net, open(colorfile, 'r') as mr:
        load_net = json.load(load_net)
        if org_file:
            #org_list = open(org_file,'r').read().strip().split('\n')
            org_list = []
            with open(org_file,'r') as fr:
                for line in fr:
                    org_list.append(line.split('\t')[3][-5:])
        else:
            org_list = load_net.keys()

        for ko in load_net:
            if org_file:
                if ko[-5:] not in set(org_list):  #
                    continue
            if ko[-5:] in primary_list:
                continue
            mr.seek(0, 0)
            node2color = dict()
            if len(load_net[ko]) == 0:
                continue
            net = load_net[ko]
            net_ = sum(net, [])
            nodes = [x.split('cpd:')[1] for x in net_ if 'cpd:' in x]
            n_map = []
            for line in mr:
                tmp_ = line.strip().split('\t')
                if len(tmp_) > 2:
                    n = tmp_[1]
                    col = tmp_[2]
                    if ';' in n:
                        #print(n)
                        n_s = n.split(';')
                        for n_1 in n_s:
                            if n_1 in nodes:
                                node2color[n_1] = col
                                n_map.append(n_1)

                    else:
                        if n in nodes:
                            node2color[n] = col
                            n_map.append(n)
                else:
                    n = tmp_[1]
                    if ';' in n:
                        #print(n)
                        n_s = n.split(';')
                        for n_1 in n_s:
                            if n_1 in nodes:
                                node2color[n_1] = "red"
                                n_map.append(n_1)
                    else:
                        if n in nodes:
                            node2color[n] = "red"
                            n_map.append(n)
            if not n_map:
                continue
            st_list, node_list = analyse_tree(net)
            draw_cpd_dag(st_list, node_list, node2color, '%s/%s.network' % (outpath, ko))
