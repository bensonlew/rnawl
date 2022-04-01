# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re,os

def change_faa(faa, dir, seqids, type):
    for seqid in seqids:
        seq_list = []
        n = 1
        for seq_record in SeqIO.parse(faa, "fasta"):
            m = re.search(r'(.*)\s+#\s+(.*)\s+#\s+(.*)\s+#\s+(.*)\s+#\s+(.*)', seq_record.description)
            if m.group(2) == seqid:
                id = ''
                str = "UN_{0:0>6}".format(n)
                if m.group(5) in ["1"]:
                    des = m.group(3) + ".." + m.group(4)
                    id = "ref|{}|gi|{}|{}".format(str, n, des)
                elif m.group(5) in ["-1"]:
                    des = "c" + m.group(3) + ".." + m.group(4)
                    id = "ref|{}|gi|{}|{}".format(str, n, des)
                seq_record.id = id
                seq_list.append(seq_record)
                n += 1
        SeqIO.write(seq_list, dir + "/" + seqid + "/" + seqid + "." + type, "fasta")


def chang_gff(gff, dir, seqids):
    for seqid in seqids:
        with open(gff, "r") as f, open(dir + "/" + seqid + "/" + seqid + ".ptt", "w") as d:
            d.write("#island_input file\n#proteins\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
            lines = f.readlines()
            n = 1
            for line in lines:
                if line.startswith(("#")):
                    continue
                else:
                    lin = line.strip().split("\t")
                    if lin[0] == seqid:
                        temp = lin[8].split(";")
                        id = re.search('ID=(.*)', temp[0])
                        len = abs((int(lin[4]) - int(lin[3]) + 1) / 3 + 1)
                        dd = [lin[3] + ".." + lin[4], lin[6], str(len), str(n), "-", id.group(1), lin[0], "-", "-"]
                        d.write("{}\n".format("\t".join(dd)))
                        n += 1



def get_seqids(genome, dir):
    seqids = []
    for seq_record in SeqIO.parse(genome, "fasta"):
        seqids.append(seq_record.id)
        os.mkdir(dir + "/" + seq_record.id)
    return seqids

def main():
    parser = OptionParser()
    parser.add_option('--s', dest='genome', metavar='[genome or metagenome file]')
    parser.add_option('--f', dest ='faa',metavar='[gene protein file]')
    parser.add_option('--n', dest ='ffn',metavar='[gene nucl file]')
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option("--o", dest="output", metavar="[out dir]")
    (options,args) = parser.parse_args()
    if not options.faa or not options.ffn or not options.gff or not options.output:
        print "python Is_gene_dir.py --s genome --f faa --n ffn --g gff --o dir"
        return
    list1 = get_seqids(options.genome, options.output)
    chang_gff(options.gff, options.output, list1)
    change_faa(options.faa, options.output, list1, "faa")
    change_faa(options.ffn, options.output, list1, "ffn")


if __name__=='__main__':
    main()