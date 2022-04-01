#-*- coding: utf-8 -*
# __author__: qingchen.zhang
from Bio import SeqIO
import os
import re
import sys
from collections import defaultdict

stat_file = sys.argv[1]
chromosome = sys.argv[2]
sample = sys.argv[3]
summary_file = sys.argv[4]
out_stat_file = sys.argv[5]
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

stat_number = 0
type_list = []
sequence_aa_dict = defaultdict(list)
clan_seq_dict = {}
complete_seq_dict = {}
in0_seq_dict = {}
stat2_path = sample + ".sample.xls"
with open(summary_file, 'r') as stat, open(stat2_path, 'w') as sw:
    stat.readline()
    sw.write("Sample\tIntegron_number\n")
    for line in stat:
        spline = line.strip().split("\t")
        if "ID_replicon" in spline:
            type_list = spline[1:]
        else:
            seq_id = spline[0]
            if int(spline[1]) != 0:
                stat_number += 1
                if seq_id not in clan_seq_dict:
                    clan_seq_dict[seq_id] = int(spline[1])
                else:
                    number = clan_seq_dict[seq_id]
                    clan_seq_dict[seq_id] = number + int(spline[1])
            if int(spline[2]) != 0:
                stat_number += 1
                if seq_id not in complete_seq_dict:
                    complete_seq_dict[seq_id] = int(spline[1])
                else:
                    number = complete_seq_dict[seq_id]
                    complete_seq_dict[seq_id] = number + int(spline[1])
            if int(spline[3]) != 0:
                stat_number += 1
                if seq_id not in in0_seq_dict:
                    in0_seq_dict[seq_id] = int(spline[1])
                else:
                    number = in0_seq_dict[seq_id]
                    in0_seq_dict[seq_id] = number + int(spline[1])
    sw.write("{}\t{}\n".format(sample, str(stat_number)))



sample_path = sample + "_sequence.fna"
stat_path = sample + ".stat.xls"
seq_path = sample + ".integron.fna"
with open(stat_file, 'r') as f, open(sample_path, 'w') as w, open(stat_path, 'w') as fw, open(seq_path, 'w') as integron:
    fw.write("Integron_ID\tLocation\tSample\tStrand\tStart\tEnd\tLength(bp)\tType\tCDS\n")
    integron_num = []
    integron_dict = {}
    cds_num = 0
    for line in f:
        if line.startswith("#"):
            continue
        elif "ID_integron" in line:
            continue
        else:
            spline = line.strip().split("\t")
            integron_id_list = spline[0].split("_")
            integron_id = integron_id_list[0] + integron_id_list[1]
            location_id = spline[1]
            sequence_id = spline[2]
            strand = spline[5]
            if strand in [1, '1']:
                start = int(spline[3])
                end = int(spline[4])
                if start < end:
                    sequence = str(seq_dict[location_id])[start-1:end]
                else:
                    sequence = str(seq_dict[location_id])[end-1:start]
            else:
                start = int(spline[4])
                end = int(spline[3])
                if start < end:
                    sequence = str(reverse_dict[location_id])[start-1:end]
                else:
                    sequence = str(reverse_dict[location_id])[end-1:start]
            w.write(">{}\n{}\n".format(sequence_id, sequence))
            if integron_id not in integron_num:
                cds_num = 0
                integron_num.append(integron_id)
                if spline[8] in ['protein']:
                    cds_num += 1
                if int(spline[3]) < int(spline[4]):
                    aa_start = int(spline[3])
                    end = int(spline[4])
                else:
                    aa_start = int(spline[4])
                    end =  int(spline[3])
            else:
                if int(spline[3]) < int(spline[4]):
                    end = int(spline[4])
                else:
                    end = int(spline[3])
                if spline[8] in ['protein']:
                    cds_num += 1
            new_integron = "In" +  integron_id_list[1]
            if strand in [1, '1']:
                new_strand = "+"
            else:
                new_strand = "-"

            ## 循环更新end length和cds_num
            integron_type = spline[10]
            if new_integron not in integron_dict:
                length = abs(end - aa_start)
                integron_list = [new_integron, location_id, sample, new_strand, aa_start, end, length, integron_type, cds_num]
                integron_dict[new_integron] = "\t".join([str(x) for x in integron_list])
            else:
                integron_list = integron_dict[new_integron].split("\t")
                start = int(integron_list[4])
                length = abs(end - start)
                integron_list[6] = length
                integron_list[5] = end
                integron_list[-1] = cds_num
                integron_dict[new_integron] = "\t".join([str(x) for x in integron_list])
    for key in integron_dict.keys():
        value = integron_dict[key]
        fw.write(value + "\n")
        value_list = value.split("\t")
        start = int(value_list[4])
        end = int(value_list[5])
        strand = value_list[3]
        location = value_list[1]
        seq_id = value_list[0]
        seq_string = ""
        if strand in ['+']:
            if location in seq_dict:
                seq_string = str(seq_dict[location])[start-1:end]
        else:
            if location in reverse_dict:
                seq_string = str(reverse_dict[location])[start-1:end]
        integron.write(">{}\n{}\n".format(seq_id, seq_string))





