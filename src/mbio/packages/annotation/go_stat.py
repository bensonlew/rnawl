# -*- coding: utf-8 -*-
# __author__ = 'qiuping'


def get_gene_go(self, go_result, gene_list, outpath):
    """
    将go_annotation注释的结果文件筛选只包含基因的结果信息,保留含有基因序列的行
    go_result:go_annotation tool运行得到的blast2go.annot或query_gos.list结果文件；
    gene_list: 只包含基因序列名字的列表
    return: gene_cog_list.xls
    """
    with open(go_result, 'rb') as c, open(outpath, 'wb') as w:
        for line in c:
            name = line.strip('\n').split('\t')[0]
            if name in gene_list:
                w.write(line)
