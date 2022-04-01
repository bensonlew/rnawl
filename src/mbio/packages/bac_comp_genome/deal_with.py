# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/29'
import os
import re
import pandas as pd
from functools import reduce

class deal(object):
    def __init__(self):
        self.order_list = ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']

    def deal_data(self, input, output):
        dataset = pd.read_table(input, sep='\t', header=0, low_memory=False)
        data = dataset.set_index('Assembly_accession')
        genome_accession = data['Genome_accession']
        gc = data['CG_content'].apply(str)
        print data['Genome_length'].head()
        genome_length = data['Genome_length'].astype(str)
        trna = data['tRNA_num'].apply(str)
        ncrna = data['ncRNA_num'].apply(str)
        rrna = data['rRNA_num'].apply(str)
        tmrna = data['tmRNA_num'].apply(str)
        gene = data['Gene_num'].apply(str)
        total_length = data['Total_length'].astype(int)
        n_length = data['N_length'].astype(int)
        special_num = data['Special_num'].astype(int)
        gene_tag = data['Gene_tag'].apply(str)
        chromo_plas = data['Plasmid/Chromosome'].apply(str)
        aa_list = [genome_accession, gc,trna,ncrna,rrna,tmrna,gene, gene_tag, chromo_plas, genome_length]
        bb_list = [total_length, n_length, special_num]
        total_list = []
        for file in aa_list:
            file1 = file.groupby(file.index).apply(','.join)
            file2 = file1.reset_index()
            total_list.append(file2)
        for num in bb_list:
            num1 = num.groupby(num.index).sum()
            num2 = num1.reset_index()
            total_list.append(num2)
        result = reduce(lambda left,right: pd.merge(left,right,on='Assembly_accession'), total_list)
        result.to_csv(output, sep='\t',index=0)

if __name__ == '__main__':
    m = deal()
    #m.deal_data('sort_head', "six_toatl_assembly_accession.xls")
    m.deal_data('sort_head', "seven_toatl_assembly_accession.xls")
