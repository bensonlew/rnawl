from Bio import SeqIO
import os
import re
import sys

stat_file = sys.argv[1]
chromosome = sys.argv[2]
sample = sys.argv[3]
gene_fna = sys.argv[4]
gene_ffn = sys.argv[5]
seq_dict = {}
reverse_dict = {}
for seq_record in SeqIO.parse(chromosome, 'fasta'):
    id = seq_record.id
    seq = str(seq_record.seq)
    reverse = str(seq_record.seq.reverse_complement())
    if id not in seq_dict:
        seq_dict[id] = seq
    if id not in reverse_dict:
        reverse_dict[id] = reverse

gene_ffn_dict = {}
for record in SeqIO.parse(gene_ffn, 'fasta'):
    id = record.id
    description = record.description
    origin_id = description.split(" ")[1]
    if id not in gene_ffn_dict:
        gene_ffn_dict[id] = origin_id

out_gene_name = sample + ".gene.fna"
with open(out_gene_name, 'w') as w:
    for record in SeqIO.parse(gene_fna, 'fasta'):
        id = record.id
        seq = str(record.seq)
        if id in gene_ffn_dict:
            w.write(">{}\t{}\n{}\n".format(id, gene_ffn_dict[id], seq))

sample_path = sample + "_sequence.fna"
sample_path2 = sample + ".stat.xls"
with open(stat_file, 'r') as f, open(sample_path, 'w') as w, open(sample_path2, 'w') as w2:
    is_num = 0
    for line in f:
        if line.startswith("#"):
            continue
        elif "seqID" in line:
            continue
        else:
            is_num += 1
            spline = re.split(r"\s+", line)
            location_id = spline[0].strip()
            is_id = "IS" + str(is_num).zfill(2)
            strand = spline[17]
            print(strand)
            is_start = int(spline[3])
            is_end = int(spline[4])
            start1 = int(spline[7])
            end1 = int(spline[8])
            start2 = int(spline[9])
            end2 = int(spline[10])
            gene_start = int(spline[15])
            gene_end = int(spline[16])
            is_left_id = is_id + "_left"
            is_right_id = is_id + "_right"
            is_gene_id = is_id + "_gene"
            is_old_id = str(location_id) + "_" + str(is_start) +"_" + str(is_end)+"_" + str(strand)
            is_old_left_id = str(location_id) + "_" + str(start1) +"_" + str(end1)+"_" + str(strand)
            is_old_rigtht_id = str(location_id) + "_" + str(start2) +"_" + str(end2)+"_" + str(strand)
            is_old_gene_id = str(location_id) + "_" + str(gene_start) +"_" + str(gene_end)+"_" + str(strand)
            if strand in ['+']:
                is_sequence = str(seq_dict[location_id])[is_start-1:is_end-1]
                is_left = str(seq_dict[location_id])[start1-1:end1-1]
                is_rigtht = str(seq_dict[location_id])[start2-1:end2-1]
                is_gene = str(seq_dict[location_id])[gene_start-1:gene_end-1]
            else:
                is_sequence = str(reverse_dict[location_id])[is_start-1:is_end-1]
                is_left = str(reverse_dict[location_id])[start1-1:end1-1]
                is_rigtht = str(reverse_dict[location_id])[start2-1:end2-1]
                is_gene = str(reverse_dict[location_id])[gene_start-1:gene_end-1]
            w.write(">{}\t{}\n{}\n".format(is_id, is_old_id, is_sequence))
            w.write(">{}\t{}\n{}\n".format(is_left_id,is_old_left_id, is_left))
            w.write(">{}\t{}\n{}\n".format(is_right_id,is_old_rigtht_id, is_rigtht))
            if is_old_gene_id in gene_ffn_dict:
                w.write(">{}\t{}\t{}\n{}\n".format(is_gene_id,is_old_gene_id, gene_ffn_dict[is_old_gene_id],is_gene))
            else:
                w.write(">{}\t{}\n{}\n".format(is_gene_id,is_old_gene_id, is_gene))
    w2.write("Sample Name\tIS No.\n")
    w2.write("{}\t{}\n".format(sample, str(is_num)))