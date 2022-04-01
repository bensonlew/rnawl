# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from __future__ import division
import re
import sys
from operator import itemgetter


class nr_stat(object):
    def __init__(self):
        self.gene_list = []

    def get_gene_list(self, gene_file):
        """
        传入基因序列文件，返回装有基因序列的列表
        gene_file: 基因序列文件
        """
        with open(gene_file, 'rb') as f:
            for line in f:
                self.gene_list.append(line.strip('\n'))
        return self.gene_list

    def detail_to_level(self, detail_file, out_dir, gene_list=None, trinity_mode=True):
        """
        对比对到nr库的结果文件query_taxons_detail.xls筛选只含有域界门纲目科属种的注释信息，去掉多余的部分（query_taxons.xls）；另外，当某一水平没有注释到相应的物种时，则该水平的物种名则命名为unclassified并加上它的父辈水平的名字。
        例如：某基因在种水平上没有注释到物种，则该种水平命名为：s__unclassified_g__xxx
        若传入gene_list，则筛选出只含有基因序列的gene_taxons_detail.xls，gene_taxons.xls。
        detail_file：blast比对到nr库的结果文件query_taxons_detail.xls
        gene_list:基因序列的列表，type：list
        trinity_mode用于在新生成的xml的queryID是去除结尾的_i(数字) 的
        out_dir：输出结果文件夹路径
        """
        query_taxons = out_dir + '/query_taxons.xls'
        gene_taxons = out_dir + '/gene_taxons.xls'
        gene_detail = out_dir + '/gene_taxons_detail.xls'
        order_list = ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']
        with open(detail_file, 'rb') as f, open(query_taxons, 'wb') as w, open(gene_taxons, 'wb') as gt, open(gene_detail, 'wb') as gd:
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
                for _index, _value in enumerate(order_list):
                    if not taxon_name[_value]:
                        if _index:
                            taxon_name[_value] = '{}__unclassified_{}'.format(_value, taxon_name[order_list[_index - 1]])
                        else:
                            taxon_name[_value] = '{}__unclassified'.format(_value)
                    species_name.append(taxon_name[_value])
                w.write('{}\t{}\n'.format(query_name, ';'.join(species_name)))
                if gene_list:
                    if query_name in gene_list:
                        if trinity_mode:
                            query_name = query_name.split('_i')[0]
                        gt.write('{}\t{}\n'.format(query_name, ';'.join(species_name)))
                        foo[0] = query_name
                        w_line = '\t'.join(foo)
                        gd.write(w_line + '\n')

    def nr_stat_info(self, tran_taxonfile, gene_list, outpath, level=7):
        """
        传入只含有域界门纲目科属种的转录本的注释信息，以及基因序列的列表，并且设置分类水平，即可统计在该分类水平的物种对应的转录本、基因的数目和百分比
        tran_taxonfile：只含有域界门纲目科属种的转录本的注释信息文件
        gene_list：基因序列的列表，type：list
        outpath：统计结果的文件路径
        level：分类水平，范围：[0,1,2,3,4,5,6,7],分别对应域界门纲目科属种
        """
        count_info = {}
        tran_num = 0
        gene_num = len(gene_list)
        with open(tran_taxonfile, 'rb') as f:
            for line in f:
                line = line.strip('\n').split('\t')
                taxon_info = line[2].split(';')[-1]
                query_name = line[0]
                tran_num += 1
                species_name = taxon_info.split("{")[0]
                if species_name in count_info:
                    count_info[species_name]['tran'] += 1
                    if query_name in gene_list:
                        count_info[species_name]['gene'] += 1
                else:
                    new_dict = {}
                    new_dict['tran'] = 1
                    if query_name in gene_list:
                        new_dict['gene'] = 1
                    else:
                        new_dict['gene'] = 0
                    count_info[species_name] = new_dict
        count_list = []
        for species in count_info:
            count_info[species]['species'] = species
            count_list.append(count_info[species])
        count_list_sorted = sorted(count_list, key=itemgetter('gene'), reverse=True)

        with open(outpath, 'wb') as w:
            w.write('taxon\ttranscripts\tunigene\ttranscript_percent\tunigene_percent\n')
            for species in count_list_sorted:
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(species['species'], species['tran'], species['gene'], '%0.4g' % (species['tran'] / tran_num * 100), '%0.4g' % (species['gene'] / gene_num * 100)))

        """
        with open(outpath, 'wb') as w:
            w.write('taxon\ttranscripts\tunigene\ttranscript_percent\tunigene_percent\n')
            for species in count_info:
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(species, count_info[species]['tran'], count_info[species]['gene'], '%0.4g' % (count_info[species]['tran'] / tran_num * 100), '%0.4g' % (count_info[species]['gene'] / gene_num * 100)))
        """

    def stats(self, detail_file, out_dir, gene_list):
        self.detail_to_level(detail_file=detail_file, out_dir=out_dir, gene_list=gene_list)
        tran_taxonfile = out_dir + '/query_taxons.xls'
        outpath = out_dir + '/nr_taxon_stat'
        tmp = {0: 'd', 1: 'k', 2: 'p', 3: 'c', 4: 'o', 5: 'f', 6: 'g', 7: 's'}
        for i in tmp:
            self.nr_stat_info(tran_taxonfile=tran_taxonfile, gene_list=gene_list, outpath=outpath + '_{}.xls'.format(tmp[i]), level=i)
