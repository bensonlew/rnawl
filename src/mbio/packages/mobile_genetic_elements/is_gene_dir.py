# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re

def change_faa(faa, outfa, dict):
    seq_list =[]
    for seq_record in SeqIO.parse(faa, "fasta"):
        m = re.search(r'(.*)\s+#\s+(.*)\s+#\s+(.*)\s+#\s+(.*)\s+#\s+(.*)', seq_record.description)
        if m.group(5) in ["1"]:
            strand = "+"
        elif m.group(5) in ["-1"]:
            strand = "-"
        id = m.group(2)+"_"+m.group(3)+"_"+m.group(4)+"_" + strand
        seq_record.id =id
        seq_list.append(seq_record)
    SeqIO.write(seq_list, outfa, "fasta")

def chang_gff(gff,outgff):
    dict = {}
    with open(gff, "r") as f,open(outgff, "w") as d:
        lines = f.readlines()
        for line in lines:
            if line.startswith(("#")):
                d.write(line)
            else:
                lin =line.strip().split("\t")
                temp =lin[8].split(";")
                id =re.search('ID=(.*)', temp[0])
                des=lin[0]+"_"+lin[3]+"_"+lin[4]+"_"+lin[6]
                temp[0]="ID=" + des
                lin[8] =";".join(temp)
                d.write("{}\n".format("\t".join(lin)))

def main():
    parser = OptionParser()
    parser.add_option('--s', dest='name', metavar='[genome or metagenome file name]')
    parser.add_option('--f', dest ='faa',metavar='[gene protein file]')
    parser.add_option('--n', dest ='ffn',metavar='[gene nucl file]')
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option("--o", dest="output", metavar="[out dir]")
    (options,args) = parser.parse_args()
    if not options.faa or not options.ffn or not options.gff or not options.output:
        print "python Is_gene_dir.py --f faa --n ffn --g gff --o dir"
        return
    id_dict = chang_gff(options.gff, options.output + "/" + options.name + ".gff")
    change_faa(options.faa,options.output+"/"+ options.name +".faa", id_dict)
    change_faa(options.ffn, options.output + "/" + options.name + ".ffn", id_dict)

if __name__=='__main__':
    main()