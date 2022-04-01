# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from BCBio import GFF
import urllib
import sys
import re

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
gene2id_handle = open(out_file + ".gene2enterz", 'w')
gene2name_handle = open(out_file + ".tran2name", 'w')
gene2des_handle = open(out_file + ".tran2des", 'w')
tran_id_list = []

for seq_record in rec_list:
    for gene_feature in seq_record.features:
        if gene_feature.type.lower() == "gene" or gene_feature.type.lower() == "pseudogene":
            gene_id = gene_feature.id
            gene_name = ""
            gene_dbref = ""
            gene_description = ""
            if 'Name' in gene_feature.qualifiers:
                gene_name = gene_feature.qualifiers['Name'][0]
                if '%' in gene_name:
                    gene_name = urllib.unquote(gene_name)
            if 'Dbxref' in gene_feature.qualifiers:
                for item in gene_feature.qualifiers['Dbxref']:
                    if re.search(r'GeneID:', item):
                        gene_dbref = item.split("ID:")[1]
                    if re.search(r'GeneID_', item):
                        gene_dbref = item.split("ID_")[1]

            for rna_feature in gene_feature.sub_features:
                tran_id = rna_feature.id
                if rna_feature.type.lower() in ["cds", 'exon']:
                    if gene_id in tran_id_list:
                        pass
                    else:
                        if 'product' in rna_feature.qualifiers:
                            gene_description = rna_feature.qualifiers['product'][0]
                        gene2name_handle.write("{}\t{}\n".format(gene_id, gene_name))
                        gene2des_handle.write("{}\t{}\n".format(gene_id, gene_description))
                        gene2id_handle.write("{}\t{}\t{}\n".format(gene_id, gene_id, gene_dbref))
                        tran_id_list.append(gene_id)
                else:
                    for exon_feature in rna_feature.sub_features:
                        if 'product' in exon_feature.qualifiers and exon_feature.qualifiers['product'][0] != "":
                            gene_description = exon_feature.qualifiers['product'][0]
                            break
                        else:
                            pass


                    gene2name_handle.write("{}\t{}\n".format(tran_id, gene_name))
                    gene2des_handle.write("{}\t{}\n".format(tran_id, gene_description))
                    gene2id_handle.write("{}\t{}\t{}\n".format(gene_id, tran_id, gene_dbref))
        else:
            pass
gene2id_handle.close()
gene2name_handle.close()
gene2des_handle.close()
