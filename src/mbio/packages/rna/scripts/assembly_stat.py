import fileinput
import re
import argparse
parser = argparse.ArgumentParser(description="file")
parser.add_argument("-tmapfile", "--tmapfile", help="*.transcripts.gtf.tmap", required=True)
parser.add_argument("-transcript_file", "--cufflinks_transcript.gtf", help="cufflinks_transcript.gtf", required=True)
parser.add_argument("-out_new_trans", "--output_trans.gtf", help="output_transcript.gtf ", required=True)
parser.add_argument("-out_new_genes", "--output_genes.gtf", help="output_genes.gtf ", required=True)
parser.add_argument("-out_old_trans", "--old_trans.gtf", help="old_trans.gtf ", required=True)
parser.add_argument("-out_old_genes", "--old_genes.gtf", help="old_genes.gtf ", required=True)
args = vars(parser.parse_args())

tmap_file = args["tmapfile"]
transcript_gtf = args["cufflinks_transcript.gtf"]
info_file = args["output_trans.gtf"]
gene_file = args["output_genes.gtf"]
old_trans = args["old_trans.gtf"]
old_genes = args["old_genes.gtf"]


def get_info_dic_from_tmap(srcfile):
    dic = dict()
    for line in fileinput.input(srcfile):
        arr = line.strip().split("\t")
        p = re.compile(r'[xijuo]')
        if re.match(p, arr[2]):
            key = arr[4]
            value = arr[2]
            dic[key] = value
    return dic


def get_gene_dic_from_tmap(srcfile):
    dic = dict()
    for line in fileinput.input(srcfile):
        arr = line.strip().split("\t")
        p = re.compile(r'[u]')
        if re.match(p, arr[2]):
            key = arr[4]
            value = arr[2]
            dic[key] = value
    return dic

candidateDic = get_info_dic_from_tmap(tmap_file)
geneDic = get_gene_dic_from_tmap(tmap_file)

p = re.compile(r'transcript_id')
fw = open(info_file, "w+")
w = open(gene_file, "w+")
old_trans = open(old_trans, "w+")
old_genes = open(old_genes, "w+")
for line in fileinput.input(transcript_gtf):
    m = re.match("#.*", line)
    if not m:
        line1 = line.strip()
        arr = line1.split("\t")
        description = arr[8]
        desc_array = description.strip().split(";")
        for tag in desc_array:
            tagpair = tag.strip().split("\"")
            if len(tagpair) == 3:
                if re.match(p, tagpair[0]):
                    if candidateDic. has_key(tagpair[1]):
                        newline1 = line1+" class_code \""+candidateDic[tagpair[1]]+"\";\n"
                        fw.write(newline1)
                    else:
                        old_trans.write(line)
                    if geneDic. has_key(tagpair[1]):
                        newline2 = line1 + " class_code \"" + geneDic[tagpair[1]] + "\";\n"
                        w.write(newline2)
                    else:
                        old_genes.write(line)


fw.close()
w.close()
old_trans.close()
old_genes.close()


