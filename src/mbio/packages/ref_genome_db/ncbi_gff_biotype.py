# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from BCBio import GFF
import urllib
import sys

in_file = sys.argv[1]
gene2biotype_file = sys.argv[2]
tran2biotype_file = sys.argv[3]

in_handle = open(in_file)
rec_list = []
for rec in GFF.parse(in_handle):
    try:
        rec_list.append(rec)
    except ValueError:
        print "坐标错误"

in_handle.close()
gene2biotype = open(gene2biotype_file, 'w')
tran2biotype = open(tran2biotype_file, 'w')

for seq_record in rec_list:
    for gene_feature in seq_record.features:
        if gene_feature.type.lower() == "gene":
            gene_id = gene_feature.id
            gene_biotype = ""
            tran_biotype = ""
            for sub_feature in gene_feature.sub_features:
                tran_id = sub_feature.id
                tran_biotype = ""
                if sub_feature.type.lower() == "mrna":
                    tran_biotype = "mRNA"
                elif sub_feature.type.lower() == 'cds':
                    tran_biotype = "mRNA"
                elif sub_feature.type.lower() == 'trna':
                    tran_biotype = "tRNA"
                elif sub_feature.type.lower() == 'rrna':
                    tran_biotype = "rRNA"
                elif sub_feature.type.lower() == 'lnc_rna':
                    tran_biotype = "lncRNA"
                elif sub_feature.type.lower() == 'antisense_rna':
                    tran_biotype = "lncRNA"
                elif sub_feature.type != "":
                    tran_biotype = "other"
                else:
                    tran_biotype = "mRNA"
                tran2biotype.write("{}\t{}\n".format(tran_id, tran_biotype))
            if 'gene_biotype' in gene_feature.qualifiers:
                gene_biotype = gene_feature.qualifiers['gene_biotype'][0]
                if gene_biotype == "protein_coding":
                    gene_biotype = "mRNA"
                elif gene_biotype == "rRNA":
                    gene_biotype = "rRNA"
                elif gene_biotype == "tRNA":
                    gene_biotype = "tRNA"
                elif gene_biotype == "lncRNA":
                    gene_biotype = "lncRNA"
                elif gene_biotype == "antisense_RNA":
                    gene_biotype = "lncRNA"
                else:
                    gene_biotype = "other"
            else:
                if tran_biotype != "":
                    gene_biotype = tran_biotype
                else:
                    gene_biotype = "mRNA"
            gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
        else:
            pass
gene2biotype.close()
tran2biotype.close()