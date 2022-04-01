# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.iofile import File
import re
from collections import defaultdict


class KeggTableFile(File):
    """
    定义kegg_table.xls格式
    """
    def __init__(self):
        super(KeggTableFile, self).__init__()
        self.gene_list = []

    def check(self):
        if super(KeggTableFile, self).check():
            return True

    def get_kegg_list(self, outdir, all_list, diff_list):
        with open(self.prop['path'], 'rb') as r, open(outdir + '/kofile', 'wb') as w, open(outdir + '/all_kofile', 'wb') as a:
            r.readline()
            head = '##ko KEGG Orthology\n##Method: BLAST Options: evalue <= 1e-05; rank <= 5\n##Summary: None\n\n#Query\tKO ID|KO name|Hyperlink\n'
            w.write(head)
            a.write(head)
            for line in r:
                line = line.strip('\n').split('\t')
                if line[0] in diff_list:
                    w.write('{}\t{}\n'.format(line[0], line[1]))
                a.write('{}\t{}|{}|{}\n'.format(line[0], line[1], line[2], line[3]))
                self.gene_list.append(line[0])
            for i in all_list:
                if i not in self.gene_list:
                    a.write('{}\tNone\n'.format(i))

    def get_query(self):
        '''
        返回基因列表 刘彬旭
        '''
        with open(self.prop['path'], 'rb') as r:
            querys = []
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                querys.append(line[0])
        return querys

    def get_pathway_koid(self):
        '''
        返回两个字典：
        ko_gene:ko id对应的geneids
        path_ko:pathway id为键，值为koid的列表；
        '''
        with open(self.prop['path'], 'rb') as r:
            ko_gene = defaultdict(list)
            path_ko = defaultdict(list)
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                ko_gene[line[1]].append(line[0])
                paths = [re.sub('path:', '', i) for i in line[-1].split(';')]
                for p in paths:
                    if p:
                        path_ko[p].append(line[1])
        return ko_gene, path_ko

    def get_gene2K(self, outdir):
        gene2K = defaultdict(lambda: [])
        with open(self.prop['path'], 'rb') as r, open(outdir + '/gene2K.info', 'wb') as w:
            r.readline()
            head = 'Query\tKo id(Gene id)\n'
            w.write(head)
            for line in r:
                if line.startswith("#"):
                    continue
                line = line.strip('\n').split('\t')
                if not line[1] in gene2K[line[0]]:
                    gene2K[line[0]].append(line[1])
            for key in gene2K.keys():
                w.write('{}\t{}\n'.format(key, ";".join(gene2K[key])))

    def get_gene2path(self, outdir):
        with open(self.prop['path'], 'rb') as r, open(outdir + '/gene2path.info', 'wb') as w:
            r.readline()
            head = 'Query\tPaths\n'
            w.write(head)
            for line in r:
                line = line.strip('\n').split('\t')
                if not line[-1]:
                    continue
                w.write('{}\t{}\n'.format(line[0], line[-1]))
