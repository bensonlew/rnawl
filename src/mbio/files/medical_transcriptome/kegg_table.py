# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last update
from biocluster.iofile import File
import re
from collections import defaultdict
from biocluster.config import Config
from mbio.packages.rna.annot_config import AnnotConfig



class KeggTableFile(File):
    """
    定义kegg_table.xls格式
    """
    def __init__(self):
        super(KeggTableFile, self).__init__()
        self.gene_list = []
        self.path2ko = Config().SOFTWARE_DIR + "/database/KEGG/pathway"

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

    def get_path2ko(self):
        '''
        获取数据库pathway ko对应编号
        '''
        path2ko_dict = dict()
        path2ko_file = self.path2ko
        with open(path2ko_file, 'r') as f:
            for line in f.readlines():
                [path, ko] = line.strip("\n").split("\t")
                path = path[-5:]
                ko = ko.split(":")[1]
                if path in path2ko_dict:
                    path2ko_dict[path].add(ko)
                else:
                    path2ko_dict[path] = set([ko])
        return path2ko_dict

    def get_pathway_koid(self, kegg_version=None):
        '''
        返回两个字典：
        ko_gene:ko id对应的geneids
        path_ko:pathway id为键，值为koid的列表；
        '''

        self.path2ko = AnnotConfig().get_file_path(db="kegg", file="pathway", version=kegg_version)

        path2ko_dict = self.get_path2ko()
        with open(self.prop['path'], 'rb') as r:
            ko_gene = defaultdict(list)
            path_ko = defaultdict(list)
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                for ko in line[1].split(";"):
                    ko_gene[ko].append(line[0])
                if len(line) >= 5:
                    paths = [re.sub('path:', '', i) for i in line[4].split(';')]
                else:
                    paths = []
                for p in paths:
                    if p:
                        kos = line[1].split(";")
                        path = p[-5:]
                        for ko in kos:
                            if path in path2ko_dict and ko in path2ko_dict[path]:
                                path_ko[p].append(ko)
        return ko_gene, path_ko

    def get_pathway_koid2(self, kegg_version=None):
        '''
        返回三个个字典：
        ko_gene:ko id对应的geneids
        path_ko:pathway id为键，值为koid的列表；
        '''
        self.path2ko = AnnotConfig().get_file_path(db="kegg", file="pathway", version=kegg_version)
        path2ko_dict = self.get_path2ko()
        with open(self.prop['path'], 'rb') as r:
            ko_gene = defaultdict(list)
            path_ko = defaultdict(list)
            gene_kegg_gene = defaultdict(list)
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                for ko in line[1].split(";"):
                    ko_gene[ko].append(line[0])
                if len(line) >= 6:
                    gene_kegg_gene[line[0]].append(line[5])
                    paths = [re.sub('path:', '', i) for i in line[4].split(';')]
                    for p in paths:
                        if p:
                            kos = line[1].split(";")
                            path = p[-5:]
                            for ko in kos:
                                if path in path2ko_dict and ko in path2ko_dict[path]:
                                    path_ko[p].append(ko)
        return ko_gene, path_ko, gene_kegg_gene


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
                if len(line) >=5 and line[4] != "":
                    w.write('{}\t{}\n'.format(line[0], line[4]))

    def get_ko2spegene(self):
        ko2species_gene = dict()
        with open(self.prop['path'], 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                if len(line) <= 5 or line[5] == "":
                    continue
                else:
                    kos = line[1].split(";")
                    spe_genes = line[5].split(";")
                    for ko in kos:
                        if ko in ko2species_gene:
                            ko2species_gene[ko] = set(list(ko2species_gene[ko]) + spe_genes)
                        else:
                            ko2species_gene[ko] = set(spe_genes)
        return ko2species_gene
