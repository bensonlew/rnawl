# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import sys
import argparse
import os
import re
import collections


class mg_taxon(object):
    def __init__(self):
        self.order_list = ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']

    def detail_to_level(self, detail_file, out_dir):
        """
        best-hit 模式，取注释结果第一条
        """
        query_list = []
        query_taxons = out_dir + '/query_taxons.xls'
        query_index = {}  # 记录每个基因比对信息来自第几个hit，0开始
        with open(detail_file, 'rb') as f, open(query_taxons, 'wb') as w:
            i = 0
            for line in f:
                i += 1
                foo = line.strip('\n').split('\t')
                query_name = foo[0]
                taxons = foo[2].split(';')
                species_name = []
                taxon_name = {'d': '', 'k': '', 'p': '', 'c': '', 'o': '', 'f': '', 'g': '', 's': ''}
                if not query_name in query_list:
                    query_list.append(query_name)
                    query_index[query_name] = 0
                    if re.search(r'OTU|ASV', query_name):
                        for taxon in taxons:
                            m = re.match(r"(.+){(.+)}$", taxon)
                            m1 = m.group(1)
                            m2 = m.group(2)
                            taxon_name = self.level_name(m1, m2, taxon_name)
                        species_name = self.add_unclassifed(taxon_name, species_name)
                        w.write('{}\t{}\n'.format(query_name, ';'.join(species_name)))
                    else:
                        new_query_n = query_name.split("_")
                        new_query_name = "_".join(new_query_n[0:len(new_query_n) - 1])
                        for taxon in taxons:
                            m = re.match(r"(.+){(.+)}$", taxon)
                            m1 = m.group(1)
                            m2 = m.group(2)
                            taxon_name = self.level_name(m1, m2, taxon_name)
                        species_name = self.add_unclassifed(taxon_name, species_name)
                        w.write('{}\t{}\n'.format(new_query_name, ';'.join(species_name)))
        return query_index

    def detail_to_level_deunclassified(self, detail_file, out_dir):
        """
        deunclassified模式，取从门到种没有unclassified的，若都有，取best_one
        """
        query_dict = collections.OrderedDict()
        query_index = {}  # 记录每个基因比对信息来自第几个hit，0开始
        tmp_index = {}  # 记录每个基因组中hit的位次
        query_taxons = out_dir + '/query_taxons.xls'
        with open(detail_file, 'rb') as f, open(query_taxons, 'wb') as w:
            i = 0
            for line in f:
                i += 1
                foo = line.strip('\n').split('\t')
                query_name = foo[0]
                taxons = foo[2].split(';')
                species_name = []
                taxon_name = {'d': '', 'k': '', 'p': '', 'c': '', 'o': '', 'f': '', 'g': '', 's': ''}
                for taxon in taxons:
                    m = re.match(r"(.+){(.+)}$", taxon)
                    m1 = m.group(1)
                    m2 = m.group(2)
                    taxon_name = self.level_name(m1, m2, taxon_name)
                species_name = self.add_unclassifed(taxon_name, species_name)
                if query_name not in tmp_index:
                    tmp_index[query_name] = 0
                else:
                    tmp_index[query_name] += 1
                if query_name not in query_dict:
                    query_dict[query_name] = ';'.join(species_name)
                    query_index[query_name] = 0
                else:
                    exis_name = query_dict[query_name]
                    query_dict[query_name] = self.judge_unclassified(exis_name, ';'.join(species_name))
                    if exis_name != query_dict[query_name]:
                        query_index[query_name] = tmp_index[query_name]
            for eachquery in query_dict:
                new_query_n = eachquery.split("_")
                new_query_name = "_".join(new_query_n[0:len(new_query_n) - 1])
                final_species_name = query_dict[eachquery]
                w.write('{}\t{}\n'.format(new_query_name, final_species_name))
        return query_index

    def detail_to_level_lca(self, detail_file, out_dir):
        """
        LCA模式,最小祖先法
        """
        query_dict = collections.OrderedDict()
        query_taxons = out_dir + '/query_taxons.xls'
        query_index = {}  # 记录每一个gene 的hit数
        with open(detail_file, 'rb') as f, open(query_taxons, 'wb') as w:
            i = 0
            for line in f:
                i += 1
                foo = line.strip('\n').split('\t')
                query_name = foo[0]
                taxons = foo[2].split(';')
                species_name = []
                taxon_name = {'d': '', 'k': '', 'p': '', 'c': '', 'o': '', 'f': '', 'g': '', 's': ''}
                for taxon in taxons:
                    m = re.match(r"(.+){(.+)}$", taxon)
                    m1 = m.group(1)
                    m2 = m.group(2)
                    taxon_name = self.level_name(m1, m2, taxon_name)
                species_name = self.add_unclassifed(taxon_name, species_name)
                if not query_dict.has_key(query_name):
                    query_dict[query_name] = species_name
                    query_index[query_name] = 1
                else:
                    query_index[query_name] += 1
                    exis_name = query_dict[query_name]
                    query_dict[query_name] = self.lca_level(exis_name, species_name)
            for eachquery in query_dict:
                final_species_name = query_dict[eachquery]
                new_query_n = eachquery.split("_")
                new_query_name = "_".join(new_query_n[0:len(new_query_n) - 1])
                if final_species_name:
                    final_species_name_str = self.lca_polishing(final_species_name)
                    w.write('{}\t{}\n'.format(new_query_name, final_species_name_str))
        return query_index

    def level_name(self, m1, m2, taxon_name):
        if m2 == 'superkingdom':
            taxon_name['d'] = 'd__' + m1
        if m2 == 'kingdom':
            taxon_name['k'] = 'k__' + m1
        if m2 == 'phylum':
            taxon_name['p'] = 'p__' + m1
        if m2 == 'class':
            taxon_name['c'] = 'c__' + m1
        if m2 == 'order':
            taxon_name['o'] = 'o__' + m1
        if m2 == 'family':
            taxon_name['f'] = 'f__' + m1
        if m2 == 'genus':
            taxon_name['g'] = 'g__' + m1
        if m2 == 'species':
            taxon_name['s'] = 's__' + m1
        return taxon_name

    def add_unclassifed(self, taxon_name, species_name):
        for _index, _value in enumerate(self.order_list):
            if not taxon_name[_value]:
                if _index > 0:
                    for j in range(_index, 0, -1):
                        if j - 1 > 0:
                            if not "unclassified" in taxon_name[self.order_list[j - 1]]:
                                taxon_name[_value] = '{}__unclassified_{}'.format(_value,
                                                                                  taxon_name[self.order_list[j - 1]])
                                break
                        else:
                            taxon_name[_value] = '{}__unclassified_{}'.format(_value,
                                                                              taxon_name[self.order_list[j - 1]])
                else:
                    taxon_name[_value] = '{}__unclassified'.format(_value)
            species_name.append(taxon_name[_value])
        return species_name

    def judge_unclassified(self, name, species_name):
        """
        字典里已有的query对应物种name中从门到种是否具有unclassified,
        若没有，则不替换，
        若有，则判断新的物种名中是否具有，如果没有，替换旧物种名，如有，则和旧物种名比保留unclassified少的
        """
        p_s = name.split(";p__")[1]
        new_p_s = species_name.split(";p__")[1]
        p_s_unclass = p_s.count("unclassified")
        new_p_s_unclass = new_p_s.count("unclassified")
        if p_s_unclass > new_p_s_unclass:
            query_species = species_name
        else:
            query_species = name
        '''
        if not "unclassified" in p_s:
            query_species = name
        else:
            if not "unclassified" in new_p_s:
                query_species = species_name
            else:
                query_species = name
        '''
        return query_species

    def lca_level(self, name, species_name):
        lca_species = []
        name_len = len(name)
        new_len = len(species_name)
        min_len = min(name_len, new_len)
        for i in range(0, min_len):
            if name[i] == species_name[i]:
                lca_species.append(name[i])
            else:
                break
        return lca_species

    def lca_polishing(self, name_list):
        """
        lca法补充齐名称，和其他两种方法保持相同格式，用于后续分析
        """
        if len(name_list) < 8:
            all_polishing = "d__inconsistent;k__inconsistent;p__inconsistent;c__inconsistent;o__inconsistent;f__inconsistent;g__inconsistent;s__inconsistent"
            name_len = len(name_list)
            others = all_polishing.split(";")[name_len:8]
            new_list = name_list + others
            final_name = ";".join(new_list)
        else:
            final_name = ";".join(name_list)
        return final_name
