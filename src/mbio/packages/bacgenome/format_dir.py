#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from optparse import OptionParser
from Bio import SeqIO
import re,os
from collections import defaultdict
import shutil

## 此脚本用于format文件，输入gff文件、faa序列文件和scaffold的序列文件

def change_faa(faa, genome, dict, list,dir):
    for i in list:
        if os.path.exists(dir +"/tmp_" + i):
            shutil.rmtree(dir +"/tmp_" + i)
        os.mkdir(dir +"/tmp_" + i)
        seq_list = []
        outfa = dir +"/tmp_" + i + "/" + i + ".prt"
        out2 = dir +"/tmp_" + i + "/" + i + ".fst"
        for seq_record in SeqIO.parse(genome, "fasta"):
            if seq_record.id == i:
                SeqIO.write(seq_record, out2, "fasta")
        id_list = []
        for tuple_j in dict[i]:
            id = tuple_j[0]
            if id not in id_list:
                id_list.append(id)
        for seq_record in SeqIO.parse(faa, "fasta"):
            if seq_record.id in id_list:
                # desc = dict[i][1].split("\t")
                desc = [x[1].split("\t") for x in dict[i] if x[0] == seq_record.id]
                if desc[0][2] in ['+']:
                    description = "# {} # {} # {} # {}".format(i, desc[0][0], desc[0][1], str(1))
                else:
                    description = "# {} # {} # {} # {}".format(i, desc[0][0], desc[0][1], str(-1))
                seq_record.description = description
                seq_list.append(seq_record)
        SeqIO.write(seq_list, outfa, "fasta")

def chang_gff(gff):
    """
    获取gff文件的对应关系
    :param gff:
    :return:
    """
    config_dict = {
        'Chr': "Chromosome" ,
        "Chr1" :"Chromosome1",
        "Chr2" :"Chromosome2",
        "Chr3" : "Chromosome3",
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
    with open(gff, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(("#")):
                continue
            elif "Gene ID" in line:
                continue
            else:
                lin =line.strip().split("\t")
                temp =lin[1].split("_ORF")
                id = lin[0].strip().split("_")
                try:
                    new_id = id[-1]
                except:
                    new_id = id[0]
                if temp[0] in config_dict:
                    scaffold_id = config_dict[temp[0]]
                else:
                    scaffold_id = temp[0]
                desc = lin[2] + "\t" + lin[3] + "\t" +lin[4]
                dict[scaffold_id].append((new_id, desc))
    return dict

def get_seqid(genome):
    seqids= []
    for seq_record in SeqIO.parse(genome, "fasta"):
        seqids.append(seq_record.id)
    return seqids

def main():
    parser = OptionParser()
    parser.add_option('--s', dest='name', metavar='[genome or metagenome file name]')
    parser.add_option('--f', dest ='faa',metavar='[gene protein file]')
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option("--o", dest="output", metavar="[out dir]")
    (options,args) = parser.parse_args()
    if not options.faa or not options.name or not options.gff or not options.output:
        print "python format_dir.py --s name --f faa --g gff --o dir"
        return
    id_dict = chang_gff(options.gff)
    list_id = get_seqid(options.name)
    change_faa(options.faa, options.name, id_dict, list_id, options.output)

if __name__=='__main__':
    main()