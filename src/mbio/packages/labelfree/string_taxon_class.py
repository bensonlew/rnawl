# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# 该脚本用于string库的物种分类, 输入文件为string库物种表格，第一次使用需要下载ncbi物种分类文件

# sed '1d' class_string.xls |awk -F "\t" '{if($1==33090){print "Plants\t"$8"\t"$11}else if($1==33208){print "Animals\t"$8"\t"$11}else if($1==4751){print "Fungi\t"$8"\t"$11}else{print "Other\t"$8"\t"$11}}'  > ppi_species.xls
# python ~/sanger_dev3/src/mbio/packages/labelfree/string_taxon_class.py species.v10.5.txt
import sys
import csv
from ete3 import NCBITaxa

ncbi = NCBITaxa()

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

def main(taxids, desired_ranks, path):
    with open(path, 'w') as csvfile:
        fieldnames = ['{}_id'.format(rank) for rank in desired_ranks]
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for taxid in taxids:
            try:
                writer.writerow(get_desired_ranks(taxid, desired_ranks))
            except Exception,e:
                unknown_dict = {'{}_id'.format(rank): 'unfind' for rank in desired_ranks}
                writer.writerow(unknown_dict)

if __name__ == '__main__':
    species_file = sys.argv[1]
    with open(species_file, 'r') as f:
        f.readline()
        taxids = [line.strip().split("\t")[0] for line in  f.readlines()]
    desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    path = 'taxids.csv'
    main(taxids, desired_ranks, path)
