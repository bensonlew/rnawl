from Bio import SeqIO
import os
import re
import sys

stat_file = sys.argv[1]
chromosome = sys.argv[2]
sample = sys.argv[3]
# gene_ffn = sys.argv[4]
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

sample_path = sample + "_sequence.fna"
sample_path2 = sample + ".stat.xls"
with open(stat_file, 'r') as f, open(sample_path, 'w') as w, open(sample_path2, 'w') as w2:
    is_num = 0
    f.readline()
    for line in f:
        if line.startswith("#"):
            continue
        elif "seqID" in line:
            continue
        else:
            is_num += 1
            spline = re.split(r"\t", line)
            location_id = spline[1].strip()
            is_id = location_id + "_" +"IS" + str(spline[0].strip().split("_")[-1]).zfill(3)
            strand = spline[4]
            is_start = int(spline[2])
            is_end = int(spline[3])
            if spline[7].strip() != "-":
                left_start = int(spline[7].strip().split("..")[0])
                left_end = int(spline[7].strip().split("..")[1])
            if spline[8].strip() != "-":
                right_start = int(spline[8].strip().split("..")[0])
                rigth_end = int(spline[8].strip().split("..")[1])
            is_left_id = is_id + "_left"
            is_right_id = is_id + "_right"
            is_gene_id = is_id + "_gene"
            if strand in ['+']:
                is_sequence = str(seq_dict[location_id])[is_start-1:is_end]
                if spline[7].strip() != "-":
                    is_left = str(seq_dict[location_id])[left_start-1:left_end]
                if spline[8].strip() != "-":
                    is_rigtht = str(seq_dict[location_id])[right_start-1:rigth_end]
            else:
                is_sequence = str(reverse_dict[location_id])[is_start-1:is_end]
                if spline[7].strip() != "-":
                    is_left = str(reverse_dict[location_id])[left_start-1:left_end]
                if spline[8].strip() != "-":
                    is_rigtht = str(reverse_dict[location_id])[right_start-1:rigth_end]
            w.write(">{}\n{}\n".format(is_id, is_sequence))
            if spline[7].strip() != "-":
                w.write(">{}\n{}\n".format(is_left_id, is_left))
            if spline[8].strip() != "-":
                w.write(">{}\n{}\n".format(is_right_id, is_rigtht))
    w2.write("Sample Name\tIS No.\n")
    w2.write("{}\t{}\n".format(sample, str(is_num)))