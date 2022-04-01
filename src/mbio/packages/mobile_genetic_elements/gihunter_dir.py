# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(genome, gff, rnt, dir):
    for seq_record in SeqIO.parse(genome, "fasta"):
        if os.path.exists(dir + "/" + seq_record.id):
            shutil.rmtree(dir + "/" + seq_record.id)
        os.mkdir(dir + "/" + seq_record.id)
        SeqIO.write(seq_record, dir + "/" + seq_record.id + '/' + seq_record.id + ".fna", "fasta")
        with open(gff, "r") as f, open(dir + "/" + seq_record.id + '/' + seq_record.id + ".ptt", "w") as d:
            d.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
            lines = f.readlines()
            n = 1
            for line in lines:
                if line.startswith(("#")):
                    continue
                else:
                    lin = line.strip().split("\t")
                    if lin[0] == seq_record.id:
                        temp = lin[8].split(";")
                        id = re.search('ID=(.*)', temp[0])
                        len = abs((int(lin[4]) - int(lin[3]) + 1) / 3 + 1)
                        dd = [lin[3] + ".." + lin[4], lin[6], str(len), str(n), "-", id.group(1), lin[0], "-", "-"]
                        d.write("{}\n".format("\t".join(dd)))
                        n += 1
        with open(rnt, "r") as f, open(dir + "/" + seq_record.id + '/' + seq_record.id + ".rnt", "w") as d:
            d.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
            lines = f.readlines()
            for line in lines:
                if line.startswith(("#", 'Location')):
                    continue
                elif line in ['\n', '\r\n']:
                    continue
                else:
                    lin = line.strip().split("\t")
                    if lin[3] == seq_record.id:
                        d.write("{}\n".format("\t".join(lin)))

def main():
    parser = OptionParser()
    parser.add_option('--s', dest='genome', metavar='[genome or metagenome file name]')
    parser.add_option('--r', dest ='rnt',metavar='[rRNA gff file]')
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option("--o", dest="output", metavar="[out dir]")
    (options,args) = parser.parse_args()
    if not options.genome or not options.rnt or not options.gff or not options.output:
        print "python Is_gene_dir.py --s genome --r rnt --g gff --o dir"
        return
    chang_gff(options.genome, options.gff, options.rnt, options.output)


if __name__=='__main__':
    main()