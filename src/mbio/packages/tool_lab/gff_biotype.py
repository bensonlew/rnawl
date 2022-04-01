# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from BCBio import GFF
import urllib
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]
print(in_file)
print(out_file)
in_handle = open(in_file)
rec_list = []
for rec in GFF.parse(in_handle):
    try:
        rec_list.append(rec)
    except ValueError:
        print("坐标错误")

in_handle.close()
gene2biotype = open(out_file, 'w')
gene_ids = dict()
gene2biotype.write("gene_id\tgene_biotype\n")

for seq_record in rec_list:
    for gene_feature in seq_record.features:
        print(gene_feature.type)
        if gene_feature.type.lower() == "gene":
            gene_biotype = ""
            if 'locus_tag' in gene_feature.qualifiers:
                gene_id = gene_feature.qualifiers['locus_tag'][0]
            else:
                gene_id = gene_feature.id
            if gene_id.startswith("gene-"):
                gene_id = gene_id.split("gene-")[1]
            print(gene_id)
            gene_ids[gene_id] = 1
            if 'gene_biotype' in gene_feature.qualifiers:
                gene_biotype = gene_feature.qualifiers['gene_biotype'][0]
                if gene_biotype == "protein_coding":
                    gene_biotype = "mRNA"
                elif gene_biotype == "rRNA":
                    gene_biotype = "rRNA"
                elif gene_biotype == "tRNA":
                    gene_biotype = "tRNA"
                elif gene_biotype == 'lncRNA':
                    gene_biotype = 'lncRNA'
                elif 'pseudogene' in gene_biotype:
                    gene_biotype = 'pseudogene'
                elif gene_biotype == 'miRNA':
                    gene_biotype = 'miRNA'
                elif gene_biotype == 'snRNA':
                    gene_biotype = 'snRNA'
                else:
                    if gene_biotype != "":
                        gene_biotype = gene_biotype
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
                        gene_biotype = "other"
                if gene_biotype == "":
                    gene_biotype = "other"
            gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
        else:
            gene_id = gene_feature.id
            if gene_id.startswith("gene-"):
                gene_id = gene_id.split("gene-")[1]
            print(gene_id)
            if gene_id not in gene_ids:
                if gene_feature.type.lower() == "trna":
                    gene_biotype = "tRNA"
                    gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
                elif gene_feature.type.lower() == "rrna":
                    gene_biotype = "rRNA"
                    gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
                elif gene_feature.type.lower() == "mirna":
                    gene_biotype = "miRNA"
                    gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
                elif gene_feature.type.lower() == "lnc_rna":
                    gene_biotype = "lncRNA"
                    gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
                elif gene_feature.type.lower() == "antisense_rna":
                    gene_biotype = "lncRNA"
                    gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
                elif gene_feature.type.lower() == "pseudogene":
                    gene_biotype = "pseudogene"
                    gene2biotype.write("{}\t{}\n".format(gene_id, gene_biotype))
                else:
                    continue
            else:
                continue
gene2biotype.close()