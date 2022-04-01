# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modified: liubinxu 20180806

import sys
import os
import re
from biocluster.config import Config

'''
import MySQLdb
import MySQLdb.cursors

db = MySQLdb.connect(host=sys.argv[2], user=sys.argv[
                     3], passwd=sys.argv[4], db="b2gdb")
cur = db.cursor(MySQLdb.cursors.DictCursor)
'''

get_relationship = ['part_of',
                    'negatively_regulates',
                    'positively_regulates',
                    'regulates',
                    'occurs_in',
                    'capable_of',
                    'capable_of_part_of', ]  # 所有的关系类型，当前仅关注这些类型
relationship_dict = {'is_a': ('#130c0e', 'solid'),
                     'part_of': ('#2a5caa', 'solid'),
                     'negatively_regulates': ('#d71345', 'solid'),
                     'positively_regulates': ('#1d953f', 'solid'),
                     'regulates': ('#ffc20e', 'solid'),
                     'occurs_in': ('#008792', 'solid'),
                     'capable_of': ('#33a3dc', 'dashed'),
                     'capable_of_part_of': ('#f36c21', 'dashed')}

class Term(object):
    """解析一个term来自obo文件的字符串"""
    def __init__(self, value):
        '''
        初始化一个go term对象
        '''
        if value:
            self.children = []
            self.levels = set()
            self.genes = []
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
        elif value_sp[0] == 'def':
            if hasattr(self, 'def'):
                self.defi.append(value_sp)
            else:
                self.defi = value_sp[1]
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
    def __init__(self, obo_fp=Config().SOFTWARE_DIR + "/database/Annotation/other2019/go.obo"):
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
        self.add_children_term()
        self.get_levels_term()

    def add_children_term(self):
        for term in self.values():
            if hasattr(term, 'is_a'):
                for is_a in term.is_a:
                    self[is_a[0].id].children.append(term.id)

    def get_levels_term(self):
        """自上而下判断所有GO的level"""
        for term in ['GO:0003674', 'GO:0005575', 'GO:0008150']:
            self[term].levels.add(1)
            for child_term in self[term].children:
                self.get_children_level(child_term, 1)

    def get_children_level(self, term, level):
        """给子节点添加所在的level"""
        level = level + 1
        self[term].levels.add(level)
        if len(self[term].children) == 0:
            return
        else:
            for children in self[term].children:
                self.get_children_level(children, level)

    def get_all_parent_id(self, term):
        """获取所有父节点的 id"""
        self.parent_list = [term]
        self.get_parent(term)
        return self.parent_list

    def get_parent(self, term):
        """获取父节点"""
        if hasattr(self[term], 'is_a'):
            parents = [is_a[0].id for is_a in self[term].is_a]
            self.parent_list.extend(parents)
            for parent in parents:
                self.get_parent(parent)
        else:
            return

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

if __name__ == '__main__':
    if len(sys.argv) >= 4:
        GO_db = Terms(obo_fp=sys.argv[3])
    else:
        GO_db = Terms()
    f = open(sys.argv[1]).read().split('\n')
    out_dir = sys.argv[2]
    d = {}
    allgo = set()
    allgene = []
    gene2go = {}
    # 获取gene_id go对应关系
    for linerecord in f:
        if len(linerecord) != 0:
            info = linerecord.split('\t')
            gene_id = info[0]
            gos = info[1].split(';')
            d[gene_id] = gos
            terms= []
            for item in gos:
                if item in GO_db:
                    terms.extend(GO_db.get_all_parent_id(item))
            gene2go[gene_id] = list(set(terms))
            for go_term in list(set(terms)):
                if go_term in GO_db:
                    GO_db[go_term].genes.append(gene_id)
            if gene_id not in allgene:
                allgene.append(gene_id)

    go1234file = open(out_dir + '/level4.stat.tsv', 'w')
    go123file = open(out_dir + '/level3.stat.tsv', 'w')
    go12file = open(out_dir + '/level2.stat.tsv', 'w')
    go1234file.write(
        'GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tGO Term (Lev4)\tGO ID (Lev4)\tSeq Number\tPercent\tSeq List\n')
    go123file.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tSeq Number\tPercent\tSeq List\n')
    go12file.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tSeq Number\tPercent\tSeq List\n')
    # x=4
    go_detail = open(out_dir + '/go_detail.xls', 'w')
    go_detail.write("gene_id\tgos\n")

    godetailfile = open(os.path.join(out_dir, 'go_level_detail.xls'), 'w')
    godetailfile.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tGO Term (Lev4)\tGO ID (Lev4)\tSeq Number\tSeq List\n')

    for go_lev1 in ['GO:0003674', 'GO:0005575', 'GO:0008150']:
        for go_lev2 in GO_db[go_lev1].children:
            if len(GO_db[go_lev2].genes) != 0:
                pcent = float(len(GO_db[go_lev2].genes)) / float(len(allgene))
                go12file.write( GO_db[go_lev1].name + '\t'
                                + GO_db[go_lev2].name + '\t' + go_lev2 + '\t'
                                + str(len(GO_db[go_lev2].genes)) + '\t' + str("%.8f" % pcent) + '\t' + ";".join(GO_db[go_lev2].genes) + '\n')

    for go_lev1 in ['GO:0003674', 'GO:0005575', 'GO:0008150']:
        for go_lev2 in GO_db[go_lev1].children:
            for go_lev3 in GO_db[go_lev2].children:
                if len(GO_db[go_lev3].genes) != 0:
                    pcent = float(len(GO_db[go_lev3].genes)) / float(len(allgene))
                    go123file.write( GO_db[go_lev1].name + '\t'
                                    + GO_db[go_lev2].name + '\t' + go_lev2 + '\t'
                                    + GO_db[go_lev3].name + '\t' + go_lev3 + '\t'
                                    + str(len(GO_db[go_lev3].genes)) + '\t' + str("%.8f" % pcent) + '\t'+ ";".join(GO_db[go_lev3].genes) + '\n')

    for go_lev1 in ['GO:0003674', 'GO:0005575', 'GO:0008150']:
        for go_lev2 in GO_db[go_lev1].children:
            for go_lev3 in GO_db[go_lev2].children:
                for go_lev4 in GO_db[go_lev3].children:
                    if len(GO_db[go_lev4].genes) != 0:
                        pcent = float(len(GO_db[go_lev4].genes)) / float(len(allgene))
                        go1234file.write( GO_db[go_lev1].name + '\t'
                                          + GO_db[go_lev2].name + '\t' + go_lev2 + '\t'
                                          + GO_db[go_lev3].name + '\t' + go_lev3 + '\t'
                                          + GO_db[go_lev4].name + '\t' + go_lev4 + '\t'
                                          + str(len(GO_db[go_lev4].genes)) + '\t' + str("%.8f" % pcent) + '\t'+ ";".join(GO_db[go_lev4].genes) + '\n')
                        godetailfile.write(GO_db[go_lev1].name + '\t'
                                           + GO_db[go_lev2].name + '\t' + go_lev2 + '\t'
                                           + GO_db[go_lev3].name + '\t' + go_lev3 + '\t'
                                           + GO_db[go_lev4].name + '\t' + go_lev4 + '\t'
                                           + str(len(GO_db[go_lev4].genes)) + '\t' + ";".join(GO_db[go_lev4].genes) + '\n')


    for gene_id,gos in gene2go.items():
        gos_show = [go_term + "(" + ",".join([str(lev) for lev in GO_db[go_term].levels]) + ")" for go_term in gos if go_term in GO_db]
        go_detail.write(gene_id + "\t" + ";".join(gos_show) + '\n')

    go1234file.close
    go123file.close
    go12file.close
    go_detail.close
