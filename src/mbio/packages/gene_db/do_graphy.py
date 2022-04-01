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
import textwrap
import xml.etree.ElementTree as ET
import sys

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
        elif value_sp[0] == 'xref':
            try:
                db, db_id = value_sp[1].split(":")
                # print db,db_id
                if hasattr(self, db):
                    getattr(self, db).append(db_id)
                else:
                    setattr(self, db, [db_id])
            except Exception as e:
                print e


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




def term_2table(xrefs=['ICD10CM','ICD9CM','NCI','MESH','OMIM','ORDO','GARD', 'EFO', 'UMLS_CUI', 'SNOMEDCT_US_2020_03_01', 'SNOMEDCT_US_2019_09_01', 'SNOMED_CT_US_2018_03_01'], obo="/mnt/ilustre/users/sanger-dev/app/database/gene_db/HumanDiseaseOntology-2020-04-20/src/ontology/doid.obo", table='do_id.xls'):
    terms = Terms(obo_fp=obo)
    with open(table, 'w') as f:
        f.write("\t".join(['DOID'] + xrefs) + "\n")
        for do_id,term in terms.items():
            f.write(do_id)
            for xref in xrefs:
                if hasattr(term, xref):
                    f.write("\t" + ",".join(getattr(term, xref)))
                else:
                    f.write("")
            f.write("\n")


def term_add_xref(term=None, xref='ORDO', obo="/mnt/ilustre/users/sanger-dev/app/database/gene_db/HumanDiseaseOntology-2020-04-20/src/ontology/doid.obo", table=None):
    terms = Terms(obo_fp=obo)
    
    with open(term, 'r') as f_in, open(table, 'w') as f_out:
        header = f_in.readline()
        f_out.write(header.strip("\n") + "\t" + xref + "\n")
        for line in f_in:
            cols = line.strip("\n").split("\t")
            orgo_ids = list()
            do_ids = cols[1].split("|")
            for do_id in do_ids:
                if do_id in terms:
                    term = terms[do_id]
                    if hasattr(term, xref):
                        orgo_ids.append(",".join(getattr(term, xref)))
                    else:
                        orgo_ids.append("")
                else:
                    orgo_ids.append("")
            f_out.write(line.strip("\n") + "\t" + "|".join(orgo_ids) + "\n")



if __name__ == '__main__':
    # term_2table()
    term_add_xref(term=sys.argv[1], table=sys.argv[2])
