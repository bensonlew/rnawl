# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re,os
import shutil
import pandas as pd
from pandas.core.frame import DataFrame


def chang_result(fa, input, sample, length, prefix):
    with open(input, "r") as f, open(prefix+ ".stat.xls", "w") as g,open(prefix+ ".sRNA.fasta", "w") as s:
        lines = f.readlines()
        num = 0
        srna_list = []
        length_srna = 0
        seq_ids = []
        for line in lines:
            if re.search("#", line):
                continue
            else:
                num += 1
                lin = re.split(r"[ ]+", line.strip())
                if lin[3] not in seq_ids:
                    seq_ids.append(lin[3])
                len1 = abs(int(lin[10]) - int(lin[9]))+1
                length_srna += len1
                lin2 = re.split(r"(\s{2})", line.strip())
                lin3 =re.split(r"([-\"]\s)", lin2[-1])
                seq1 = get_seq(fa, lin[3], lin[9], lin[10], lin[11])
                srna_list.append([lin[3], lin[2], lin[1], lin3[2], lin[9], lin[10], len1, lin[11], lin[17], lin[16], seq1])
        ce = float(length_srna) / int(length) *100
        g.write("sRNA Num\tTotal Len (bp)\tIn Genome (%)\n")
        g.write("{}\t{}\t{}\n".format(str(num), str(length_srna), str(ce)))
        data = DataFrame(srna_list)
        data.columns = ['location', "rfam_id", "family", "family_des", "start", "end", "length", "strand", "evalue", "score", "seq"]
        list2= []
        for i in seq_ids:
            data1 = data[data["location"] == i]
            data1["start"] = data1["start"].astype(int)
            data1.sort_values(["start"], inplace=True)
            srnaids = []
            for i in range(1,data1.shape[0]+1):
                if 10 <= len(srna_list) < 100:
                    srnaid = "sRNA"+'{:d}'.format(i).zfill(2)
                elif 100 <= len(srna_list) < 1000:
                    srnaid = "sRNA"+'{:d}'.format(i).zfill(3)
                else:
                    srnaid = "sRNA"+'{:d}'.format(i).zfill(2)
                srnaids.append(srnaid)
            data1['srna_id'] = srnaids
            list2.append(data1)
        all_data = pd.concat(list2)
        all_data['sample'] = sample
        all_data.to_csv(prefix + ".sRNA.xls", sep='\t', header=True, index=False, columns=["sample",'location', "srna_id","rfam_id", "family", "family_des", "start", "end", "length", "strand", "evalue", "score", "seq"])
        for i in all_data.ix[:, ['location', "srna_id", "seq"]].values.tolist():
            s.write(">{}\n{}\n".format(i[0] + "_" + i[1], i[2]))

def get_seq(fa, loction, start, end, strand):
    seq =''
    for seq_record in SeqIO.parse(fa, "fasta"):
        if loction == seq_record.id:
                if strand == "+":
                    seq = str(seq_record.seq[int(start) - 1:int(end)])
                elif strand == "-":
                    seq1 = str(seq_record.seq[int(end) - 1:int(start)])
                    seq = dna_reverse(dna_complement(seq1))
    return seq

def dna_complement(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()

def dna_reverse(sequence):
    sequence = sequence.upper()
    return sequence[::-1]



def main():
    parser = OptionParser()
    parser.add_option('--fa', dest='fasta', metavar='[genome file]')
    parser.add_option('--in', dest='srna_tblout', metavar='[sRNA tblout file]')
    parser.add_option('--s', dest='sample', metavar='[sample name]')
    parser.add_option('--len', dest='genome_length', metavar='[genome length]')
    parser.add_option("--o", dest="prefix", metavar="[out file prefix]")
    (options,args) = parser.parse_args()
    if not options.fasta or not options.sample or not options.genome_length or not options.srna_tblout or not options.prefix:
        print "python enzyme_type.py --fa fasta --in tblout --s sample_name --len 78899 --o prefix "
        return
    chang_result(options.fasta, options.srna_tblout, options.sample, options.genome_length, options.prefix)

if __name__=='__main__':
    main()
