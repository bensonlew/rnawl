# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from BCBio import GFF
import urllib
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

in_handle = open(in_file)
rec_list = []
for rec in GFF.parse(in_handle):
    try:
        rec_list.append(rec)
    except ValueError:
        print "坐标错误"

in_handle.close()
gene2biotype = open(out_file, 'w')

for seq_record in rec_list:
    for gene_feature in seq_record.features:
        if gene_feature.type.lower() in ["gene", 'mrna']:
            gene_id = gene_feature.id
            gene_biotype = ""
            if 'locus_tag' in gene_feature.qualifiers:
                gene_id = gene_feature.qualifiers['locus_tag'][0]
            if 'gene_biotype' in gene_feature.qualifiers:
                gene_biotype = gene_feature.qualifiers['gene_biotype'][0]
                if gene_biotype == "protein_coding":
                    gene_biotype = "mRNA"
                elif gene_biotype == "rRNA":
                    gene_biotype = "rRNA"
                elif gene_biotype == "tRNA":
                    gene_biotype = "tRNA"
                else:
                    gene_biotype = "other"
            else:
                for sub_feature in gene_feature.sub_features:
                    if  sub_feature.type.lower() == 'cds' or sub_feature.type.lower() == 'mrna':
                        gene_biotype = "mRNA"
                    elif sub_feature.type.lower() == 'trna' and gene_biotype == "":
                        gene_biotype = "tRNA"
                    elif sub_feature.type.lower() == 'rrna' and gene_biotype == "":
                        gene_biotype = "rRNA"
                    elif sub_feature.type != "" and gene_biotype == "":
                        gene_biotype = "other"
                    else:
                        gene_biotype = "mRNA"
                if gene_biotype == "":
                    gene_biotype = "mRNA"
            gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
        elif gene_feature.type.lower() == "trna":
            gene_id = gene_feature.id
            gene_biotype = "tRNA"
            gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
        elif gene_feature.type.lower() == "rrna":
            gene_id = gene_feature.id
            gene_biotype = "rRNA"
            gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
        elif gene_feature.type.lower() == "mirna":
            gene_id = gene_feature.id
            gene_biotype = "miRNA"
            gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
        elif gene_feature.type.lower() == "lnc_rna":
            gene_id = gene_feature.id
            gene_biotype = "lncRNA"
            gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
        else:
            pass
gene2biotype.close()