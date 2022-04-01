# import cPickle as cpickle
import sqlite3
# nr2des_dict = {}
import sys
import csv

g2t2p = sys.argv[1]
idmapping = sys.argv[2]
match_lines = sys.argv[3]
get_lines = sys.argv[4]
out_file = sys.argv[5]


trans_list = set()
with open(g2t2p) as f:
    for line in f:
        cols = line.strip().split("\t")
        trans_list.add(cols[1])
        if cols[1].startswith("rna-"):
            trans_list.add(cols[1].split("rna-")[1])

tran2id = dict()

with open(idmapping, 'r') as idmapping_file:
    for id_dict in csv.DictReader(idmapping_file, delimiter='\t'):
        # print id_dict
        match_values = set([id_dict.get(key, "") for key in match_lines.split(",")])
        inter = list(match_values.intersection(trans_list))

        # print inter

        if len(inter) > 0:
            tran2id[inter[0]] = id_dict


with open(g2t2p, 'r') as f, open(out_file, 'w') as fw:
    fw.write("gene_id\ttranscript_id\tprotein_id\t" + get_lines.replace(",", "\t") + "\n")
    for line in f:
        cols = line.strip("\n").split("\t")
        if len(cols) == 2:
            cols.append("")
        if cols[1] in tran2id:
            id_dict = tran2id[cols[1]]
            id_values = [id_dict.get(x, "") for x in get_lines.split(",")]
            fw.write("\t".join(cols + id_values) + "\n")
        elif cols[1].split("rna-")[-1] in tran2id:
            ncbi_tran = cols[1].split("rna-")[-1]
            id_dict = tran2id[ncbi_tran]
            id_values = [id_dict.get(x, "") for x in get_lines.split(",")]
            fw.write("\t".join(cols + id_values) + "\n")
