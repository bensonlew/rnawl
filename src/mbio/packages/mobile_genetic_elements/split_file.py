# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re,os
import shutil

def change_fa(faa, list, type, dir):
    for i in list:
        seqs = []
        for seq_record in SeqIO.parse(dir + "/" + i+"/"+i+".fna", "fasta"):
            seqs.append(seq_record.id)
        seq_list = []
        for seq_record in SeqIO.parse(faa, "fasta"):
            m = re.search(r'(.*)\s+#\s+(.*)\s+#\s+(.*)\s+#\s+(.*)\s+#\s+(.*)', seq_record.description)
            if m.group(2) in seqs:
                seq_list.append(seq_record)
        outfa = dir + "/" + i+"/"+i+"." + type
        SeqIO.write(seq_list, outfa, "fasta")


def chang_gff(gff,list, num, type, dir):
    for i in list:
        seqs = []
        for seq_record in SeqIO.parse(dir + "/" + i +"/"+i+".fna", "fasta"):
            seqs.append(seq_record.id)
        outgff = dir + "/" + i+"/"+i+"." + type
        with open(gff, "r") as f, open(outgff, "w") as d:
            lines = f.readlines()
            for line in lines:
                if line.startswith(("#")):
                    continue
                elif line.startswith(("Location")):
                    d.write(line)
                else:
                    lin = line.strip().split("\t")
                    if lin[num] in seqs:
                        d.write("\t".join(lin) + "\n")


def split_genome(genome, name, dir):
    list_name = []
    seq_list =[]
    n =1
    for seq_record in SeqIO.parse(genome, "fasta"):
        if len(seq_record.seq) >= 2000:
            seq_list.append(seq_record)
        if len(seq_list) == 100:
            if os.path.exists(dir + "/" + name + "__" + str(n)):
                shutil.rmtree(dir + "/" + name + "__" + str(n))
            os.mkdir(dir + "/" + name + "__" + str(n))
            outfa = dir + "/" + name + "__" + str(n) + "/" + name + "__" + str(n) + ".fna"
            list_name.append(name + "__" + str(n))
            SeqIO.write(seq_list, outfa, "fasta")
            n +=1
            seq_list =[]
    if len(seq_list) >=1:
        if os.path.exists(dir + "/" + name + "__" + str(n)):
            shutil.rmtree(dir + "/" + name + "__" + str(n))
        os.mkdir(dir + "/" + name + "__" + str(n))
        outfa = dir + "/" + name + "__" + str(n) + "/" + name + "__" + str(n) + ".fna"
        list_name.append(name + "__" + str(n))
        SeqIO.write(seq_list, outfa, "fasta")
    return list_name

def main():
    parser = OptionParser()
    parser.add_option('--s', dest='genome', metavar='[metagenome file name]')
    parser.add_option('--a', dest='name', metavar='[sample name]')
    parser.add_option('--f', dest ='faa',metavar='[gene protein file]')
    parser.add_option('--n', dest ='ffn',metavar='[gene nucl file]')
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option('--t', dest='trna', metavar='[gene tRNA file]')
    parser.add_option("--o", dest="output", metavar="[out dir]")
    (options,args) = parser.parse_args()
    if not options.genome or not options.name or not options.faa or not options.ffn or not options.gff or not options.trna or not options.output:
        print "python split_file.py --s genome --a name --f faa --n ffn --g gff --t trna --o dir"
        return
    id_dict = split_genome(options.genome, options.name, options.output)
    change_fa(options.faa, id_dict, "faa", options.output)
    change_fa(options.ffn, id_dict, "ffn", options.output)
    chang_gff(options.gff, id_dict, 0, "ptt", options.output)
    chang_gff(options.trna, id_dict, 3, "rnt", options.output)

if __name__=='__main__':
    main()