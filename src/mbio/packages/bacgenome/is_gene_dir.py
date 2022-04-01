# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from optparse import OptionParser
from Bio import SeqIO
import re
import os
from collections import defaultdict

def get_seqid(genome):
    seqids= []
    for seq_record in SeqIO.parse(genome, "fasta"):
        seqids.append(seq_record.id)
    return seqids


def change_faa(faa, genome, dict, list,dir, genome_fa):
    outfa = open(dir +"/" + genome + ".ffn", 'w')
    out2 = open(dir +"/" + genome  + ".faa", 'w')

    for i in list:
        id_list = []
        for tuple_j in dict[i]:
            id = tuple_j[0]
            if id not in id_list:
                id_list.append(id)
        for seq_record in SeqIO.parse(faa, "fasta"):
            if seq_record.id in id_list:
                # desc = dict[i][1].split("\t")
                desc = [x[1].split("\t") for x in dict[i] if x[0] == seq_record.id]
                print(desc)
                if desc[0][2] in ['+']:
                    description = "{}_{}_{}_{} {} # {} # {} # {} # {}".format(i, desc[0][0], desc[0][1], str("+"), seq_record.id, i, desc[0][0], desc[0][1], str(1))
                else:
                    description = "{}_{}_{}_{} {} # {} # {} # {} # {}".format(i, desc[0][0], desc[0][1], str("-"), seq_record.id, i, desc[0][0], desc[0][1], str(-1))
                out2.write(">{}\n{}\n".format(description, seq_record.seq))

        for seq_record in SeqIO.parse(genome_fa, "fasta"):
            if seq_record.id in id_list:
                if desc[0][2] in ['+']:
                    description = "{}_{}_{}_{} {} # {} # {} # {} # {}".format(i, desc[0][0], desc[0][1], str("+"), seq_record.id, i, desc[0][0], desc[0][1], str(1))
                else:
                    description = "{}_{}_{}_{} {} # {} # {} # {} # {}".format(i, desc[0][0], desc[0][1], str("-"), seq_record.id, i, desc[0][0], desc[0][1], str(-1))
                outfa.write(">{}\n{}\n".format(description, seq_record.seq))



def chang_gff(gff,outgff):
    config_dict = {
        'Chr': "Chromosome" ,
        "Chr1" :"Chromosome1",
        "Chr2" :"Chromosome2",
        "Chr3" : "Chromosome3",
        'p': "Plasmid",
        'pA':"PlasmidA",
        'pB':"PlasmidB",
        'pC':"PlasmidC",
        'pD':"PlasmidD",
        'pE':"PlasmidE",
        'pF':"PlasmidF",
        'pG':"PlasmidG",
        'pH':"PlasmidH",
        'pI':"PlasmidI",
        'pJ':"PlasmidJ",
        'pK':"PlasmidK",
        'pL':"PlasmidL",
        'pM':"PlasmidM",
        'pN':"PlasmidN",
        'pO':"PlasmidO"
    }
    dict = defaultdict(list)
    with open(gff, "r") as f,open(outgff, "w") as d:
        lines = f.readlines()
        for line in lines:
            if line.startswith(("#")):
                d.write(line)
            elif line.startswith("Gene ID"):
                d.write("##gff-version 3\n")
            else:
                temp = ["0", "1"]
                lin =line.strip().split("\t")
                location = lin[1].split("_ORF")[0]
                if location in config_dict:
                    scaffold = config_dict[location]
                else:
                    scaffold = location
                gene = lin[0].strip().split("_")
                try:
                    new_gene_id = gene[-1]
                except:
                    new_gene_id = gene[0]
                orign_prefix = "CDS"
                start = lin[2]
                end = lin[3]
                strand = lin[4]
                new_line = [scaffold, new_gene_id, orign_prefix, start, end, ".", strand, '.']
                des=scaffold+"_"+lin[2]+"_"+lin[3]+"_"+lin[4]
                temp[0]="ID=" + des
                temp[1]="Name=" + new_gene_id
                prefix =";".join(temp)
                new_line.append(prefix)
                d.write("{}\n".format("\t".join(new_line)))
                desc = lin[2] + "\t" + lin[3] + "\t" +lin[4]
                dict[scaffold].append((new_gene_id, desc))
    return dict


def main():
    parser = OptionParser()
    parser.add_option('--s', dest='name', metavar='[genome or metagenome file]')
    parser.add_option('--f', dest ='faa',metavar='[gene protein file]')
    parser.add_option('--n', dest ='ffn',metavar='[gene nucl file]')
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option("--o", dest="output", metavar="[out dir]")
    parser.add_option("--a", dest="smaple", metavar="[sample name]")
    (options,args) = parser.parse_args()
    if not options.faa or not options.ffn or not options.gff or not options.output:
        print "python Is_gene_dir.py --f faa --n ffn --g gff --o dir --a sample"
        return
    id_dict = chang_gff(options.gff, options.output + "/" + options.smaple + ".gff")
    list_id = get_seqid(options.name)
    change_faa(options.faa, options.smaple, id_dict, list_id, options.output, options.ffn)


if __name__=='__main__':
    main()