# -*- coding: utf-8 -*-
"""
脚本说明:脚本基于goatools:https://github.com/tanghaibao/goatools的obo_parser.py脚本改写，用于读取Gene Ontology的obo类型文件，存储信息。
obo文件的测试来源是go-plus.obo:http://www.geneontology.org/ontology/extensions/go-plus.obo
author: shenghe
lastmodified:20160601
"""


import os
import re
import pygraphviz as pgv
from pygraphviz import *
from colour import Color
from numpy import linspace
from math import log10

get_relationship = ['part_of', 'negatively_regulates', 'positively_regulates', 'regulates', 'occurs_in', 'capable_of',
                    'capable_of_part_of', ]  # 所有的关系类型，当前仅关注这些类型
relationship_dict = {'is_a': ('#130c0e', 'solid'),
                     'part_of': ('#2a5caa', 'solid'),
                     'negatively_regulates': ('#d71345', 'solid'),
                     'positively_regulates': ('#1d953f', 'solid'),
                     'regulates': ('#ffc20e', 'solid'),
                     'occurs_in': ('#008792', 'solid'),
                     'capable_of': ('#33a3dc', 'dashed'),
                     'capable_of_part_of': ('#f36c21', 'dashed')}  # 所有关系对应的连线的颜色和类型


class Term(object):
    """解析一个term来自obo文件的字符串"""
    def __init__(self, value):
        if value:
            self.parse_term_string(value=value)
        else:
            raise Exception('term字符串不能为空')


    def parse_term_string(self, value):
        if isinstance(value, str) or isinstance(value, unicode):
            value_sp = value.strip().split('\n')
            if len(value) == 1:
                raise Exception('错误的字符串格式，不能只有一行')
            else:
                for record in value_sp:
                    self.parse_one_record(record)

    def parse_one_record(self, value):
        value_sp = re.split(': ', value.strip(), maxsplit=1)
        self.is_obsolete = False
        if value_sp[0] == 'id':
            if hasattr(self, 'id'):
                print('WARNING:存在重复ID:{} {}'.format(self.id, value_sp[1]))
            self.id = value_sp[1]
        elif value_sp[0] == 'name':
            if hasattr(self, 'name'):
                print('WARNING:存在重复name:{} /// {}'.format(self.name, value_sp[1]))
            else:
                self.name = value_sp[1]
        elif value_sp[0] == 'alt_id':
            if hasattr(self, 'alt_id'):
                self.alt_id.append(value_sp)
            else:
                self.alt_id = [value_sp[1], ]
        elif value_sp[0] == 'namespace':
            if hasattr(self, 'namespace'):
                print('WARNING:存在重复的namespace:{} {}'.format(self.namespace, value_sp))
            else:
                self.namespace = value_sp[1]
        elif value_sp[0] == 'is_a':
            self.parse_is_a(value_sp[1])
        elif value_sp[0] == 'is_obsolete' and value_sp[1] == 'true':
            self.is_obsolete = True
        elif value_sp[0] == 'relationship':
            self.parse_relationship(value_sp[1])


    def parse_is_a(self, value):
        value_sp = value.split(' ! ')
        is_a_sp = value_sp[0].split(' ')
        if len(is_a_sp) == 2:
            if is_a_sp[1] == '{is_inferred="true"}':
                if hasattr(self, 'is_a'):
                    self.is_a.append([is_a_sp[0], True])
                else:
                    self.is_a = [[is_a_sp[0], True], ]
            else:
                print('WARNING:存在不是inferred的is_a属性:{}'.format(is_a_sp[1]))
        else:
            if hasattr(self, 'is_a'):
                self.is_a.append([is_a_sp[0], False])
            else:
                self.is_a = [[is_a_sp[0], False], ]

    def parse_relationship(self, value):
        value_sp = re.split(' ! ', value)
        relation = value_sp[0].split(' ')
        if len(relation) == 3 or 2:
            if relation[0] in get_relationship:
                if hasattr(self, 'relationship'):
                    if relation[0] in self.relationship:
                        self.relationship[relation[0]].append(relation[1])
                    else:
                        self.relationship[relation[0]] = [relation[1]]
                else:
                    self.relationship = {relation[0]: [relation[1]]}
            pass
        elif len(relation) == -2:
            pass
        else:
            print('WARNING:relationship中存在格式问题:{}'.format(value))

    def get_all_parents_edge(self, relationship=True, inferred=True):
        all_parent_edges = {}
        all_orient_terms = []
        if hasattr(self, 'is_a'):
            for is_a in self.is_a:
                if inferred:
                    all_orient_terms = [(i[0], 'is_a') for i in self.is_a]
                else:
                    for is_a in self.is_a:
                        if not is_a[1]:
                            all_orient_terms.append((is_a[0], 'is_a'))
        if relationship and hasattr(self, 'relationship'):
            for key, terms in self.relationship.iteritems():
                for term in terms:
                    all_orient_terms.append((term, key))
        for p in all_orient_terms:
            if p[1] in all_parent_edges:
                all_parent_edges[p[1]].add((self, p[0]))
            else:
                all_parent_edges[p[1]] = set([(self, p[0])])
            p_edges = p[0].get_all_parents_edge()
            for key in p_edges:
                if key in all_parent_edges:
                    all_parent_edges[key] |= p_edges[key]
                else:
                    all_parent_edges[key] = p_edges[key]
        return all_parent_edges


class Terms(dict):
    """解析GO库的OBO文件，生成一个包含所有terms的对象，类似于一个可以查询的字典"""
    def __init__(self, obo_fp="C:\\Users\\sheng.he.MAJORBIO\\Desktop\\goatools-master\\go-plus.obo"):
        """初始化所有的terms"""
        self.obo_fp = obo_fp
        self.parse_obo_file()


    def parse_obo_file(self):
        if not os.path.exists(self.obo_fp):
            raise Exception('提供的OBO文件不存在:{}'.format(self.obo_fp))
        if os.path.getsize(self.obo_fp) / 1024 / 1024 > 100:
            raise Exception('文件过大，无法操作，可能传入的文件并非GO库的OBO文件')
        with open(self.obo_fp, 'r') as obo:
            line = obo.read()
            line_split = re.split('\\n\[Term\]\\n', line)
            relations = re.split('\\n\[Typedef\]\\n', line_split[-1])
            terms = line_split[1:-1]
            terms.append(relations[0])
            obo_info = line_split[0]  # obo文件的前面为obo文件的基本说明，暂时不予解析使用
            relations = relations[1:]
            for one in terms:
                t = Term(one)
                self[t.id] = t
        self.orient_object_term()

    def orient_object_term(self):
        """修改term中的is_a和relationship的关系指向为term对象，而不是ID"""
        for term in self.values():
            if hasattr(term, 'is_a'):
                for is_a in term.is_a:
                    is_a[0] = self[is_a[0]]
            if hasattr(term, 'relationship'):
                for key, orient_terms in term.relationship.iteritems():
                    term.relationship[key] = [self[i] for i in orient_terms]

    def get_term(self, term_id):
        """通过ID获取一个term对象"""
        return self[term_id]

    def get_terms(self, terms_id):
        """通过IDs列表获取terms对象列表"""
        return [self.get_term(i) for i in terms_id]

    def get_all_parents_edge(self, terms, relationship=True, inferred=True):
        """获取一些或者一个terms的所有关系指向"""
        all_parent_edges = {}
        if isinstance(terms, list):
            if isinstance(terms[0], Term):
                pass
            else:
                terms = self.get_terms(terms)
        else:
            if isinstance(terms, Term):
                terms = [terms]
            else:
                terms = [self.get_term(terms)]
        for term in terms:
            term_edges = term.get_all_parents_edge()
            for key in term_edges:
                if key in all_parent_edges:
                    all_parent_edges[key] |= term_edges[key]
                else:
                    all_parent_edges[key] = term_edges[key]
        return all_parent_edges


def get_color(values, steps=100):
    all_red_colors = list(Color.range_to(Color('yellow'), Color('red'), steps - 1))
    all_blue_colors = list(Color.range_to(Color('green'), Color('yellow'), steps - 1))
    # print("vvvvvvvvvvalusssssss")
    # print(values)
    new_values = []
    for i in values:
        if i == 0:
            new_values.append(100)
        elif i == 1.0:
            new_values.append(1)
        else:
            new_values.append(-log10(i))
    red_range = linspace(1.3, 10, steps)
    blue_range = linspace(0, 1.3, steps)

    def get_range(value, min_i, max_i, target):
        index = (min_i + max_i) / 2
        if value > target[index]:
            if max_i - index == 1:
                return index
            else:
                return get_range(value, index, max_i, target)
        elif value < target[index]:
            if index - min_i == 1:
                return min_i
            else:
                return get_range(value, min_i, index, target)
    colors = []

    for i in new_values:
        if i > 10:
            colors.append(all_red_colors[-2])
        elif i < 0:
            colors.append(all_blue_colors[1])
        elif i == 1:
            colors.append(Color('grey'))
        elif i <= 1.3:
            colors.append(all_blue_colors[get_range(i, 0, steps - 1, blue_range)])
        elif i <= 10:
            colors.append(all_red_colors[get_range(i, 0, steps - 1, red_range)])
        else:
            raise Exception('未知异常')
    print("cccccccccolorss")
    # print [i.hex for i in colors]
    return [i.hex for i in colors]

# added by gdq
def get_node_bg_colors(pvalues, pvalue_cutoff=0.05):
    from matplotlib import colors
    import math
    negative_log10pvalues = list()
    for each in pvalues:
        if each == 0 or each == 0.0:
            negative_log10pvalues.append(11)
        else:
            negative_log10pvalues.append(-math.log(each, 10))
    # color_pool = ['blue', 'red']
    color_pool = ["#FF556A", "#80FB39", "#35F3F8", "#7503EB"]
    # color_pool = list(reversed(color_pool))
    cmap = colors.LinearSegmentedColormap.from_list('my_color', color_pool, 256, )
    # color_convert = colors.Normalize(vmin=-math.log(pvalue_cutoff, 10), vmax=10, clip=False)
    color_convert = colors.Normalize(vmin=0, vmax=1, clip=True)
    # cmap.set_under('gray')
    # cmap.set_over('red')
    # color_list = cmap([color_convert(x) for x in negative_log10pvalues])
    color_list = cmap([color_convert(x) for x in pvalues])
    hex_color_list = list()
    for each in color_list:
        each[3] = 0.6
        hex_color_list.append(colors.rgb2hex(each))
    return hex_color_list


def draw_GO(GOs, out='GO_lineage', obo="/mnt/ilustre/users/sanger-dev/app/database/GO/go-plus.obo"):
    """"""
    if isinstance(GOs, list):
        # GOs = GOs[:20]
        terms_colors = [(i, '#ffce7b') for i in GOs]
    if isinstance(GOs, dict):
        temp_tuple = [[i, GOs[i]] for i in GOs]
        p_value = [i[1] for i in temp_tuple]
        # temp_colors = get_color(p_value)
        temp_colors = get_node_bg_colors(p_value)
        terms_colors = zip([i[0] for i in temp_tuple], temp_colors)
        GOs = GOs.keys()
    terms = Terms(obo_fp=obo)

    G = pgv.AGraph()
    go_relation = terms.get_all_parents_edge(terms=GOs)
    for relation in go_relation:
        for edge in go_relation[relation]:
            G.add_edge(edge[1].id + r'\n' + edge[1].name,
                       edge[0].id + r'\n' + edge[0].name,
                       label=relation, color=relationship_dict[relation][0],
                       style=relationship_dict[relation][1],
                       minlen=1.5, arrowsize=1.3, penwidth=1.5)
    G.add_nodes_from([terms[i].id + r'\n' + terms[i].name for i in GOs])
    G.graph_attr.update(dpi="180")
    G.node_attr.update(shape="box", style="rounded,filled", fillcolor="#FFFFFF")
    # G.node_attr.update(shape="box", style="rounded,filled", fillcolor="#FFFFFF")
    G.edge_attr.update(dir="back")
    for i in terms_colors:
        node = G.get_node(terms[i[0]].id + r'\n' + terms[i[0]].name)
        node.attr['fillcolor'] = i[1]
    environ = os.environ['LD_LIBRARY_PATH']
    # 重置环境， app目录下的环境会报错
    os.environ['LD_LIBRARY_PATH'] = "/usr/lib/:/usr/lib64:" + os.environ['LD_LIBRARY_PATH']
    G.draw(out + '.png', prog="dot")
    os.environ['LD_LIBRARY_PATH'] = environ




if __name__ == '__main__':
    # recs = ['GO:0046466', 'GO:0030149', 'GO:0006643', 'GO:0006665', 'GO:0016042', 'GO:0044242', 'GO:0006950',
    #         'GO:0006664', 'GO:0006687', 'GO:0044699', 'GO:0044763', 'GO:0044710', 'GO:0008150', 'GO:0006952',
    #         'GO:0006629', 'GO:0045087', 'GO:0003674', 'GO:0070776', 'GO:0070775', 'GO:0006955']
    # p_bonferroni = [5.58E-09, 5.58E-09, 3.05E-08, 3.05E-08, 3.33E-08, 6.06E-08, 7.82E-08, 8.35E-08, 8.35E-08, 1.99E-07,
    #                 2.86E-07, 3.40E-07, 3.45E-07, 4.72E-07, 7.47E-07, 8.18E-07, 9.87E-07, 1.02E-06, 1.02E-06, 2.29E-06]
    # my_test = dict(zip(recs, p_bonferroni))
    # draw_GO(my_test, obo="C:\\Users\\sheng.he.MAJORBIO\\Desktop\\goa\\go-basic.obo")
    # draw_GO(recs, out='GO_lineage_1', obo="C:\\Users\\sheng.he.MAJORBIO\\Desktop\\goa\\go-basic.obo")
    my_test = {'GO:1902305': 0.000651805904284, 'GO:0019836': 0.000704504425152,
               'GO:0044179': 0.000704504425152, 'GO:0051801': 0.000704504425152,
               'GO:0005576': 0.000724690563627, 'GO:0001897': 0.000704504425152, 'GO:0050828': 0.000704504425152,
               'GO:0052331': 0.000704504425152, 'GO:0001907': 0.000704504425152, 'GO:0044004': 0.000704504425152}
    draw_GO(my_test, out= "sj", obo= "/mnt/ilustre/users/sanger-dev/app/database/GO/go-basic.obo")
