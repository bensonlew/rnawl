# -*- coding: utf-8 -*-
from collections import defaultdict
from collections import OrderedDict
import os
import os
import re
import logging
import numpy as np
import pandas as pd
import argparse
import subprocess

def extract_pos(gene_list,gtf_path,out_file):
    id2pos = defaultdict(list)
    with open(gtf_path,"r") as gtf:
        for line in gtf.readlines():
            if not line.startswith("#"):
                line = line.strip().split("\t")
                type, start, end, score, strand, phase, attributes = line[2], line[3], line[4], line[5], line[6], line[
                    7], line[8]
                if type == "gene":
                    gene_id = re.search(r'gene_id \"(.+?)\";', attributes).group(1)
                    id2pos[gene_id] = [start, end]
    print(len(id2pos))
    with open(out_file,"w") as w,open(gene_list,"r") as g:
        for line in g.readlines():
            gene_id = line.strip()
            if gene_id in id2pos:
                w.write(gene_id+"\t"+"\t".join(id2pos[gene_id])+"\n")
            else:                       #正常情况不会出现这种情况，因为没有位置的基因star_fusion是不会输出的，这里是为了防止万一所以怼这样一波
                w.write(gene_id + "\t" + "\t".join(["unkown","unkown"])+"\n")


def extract_all_gene_pos_old(gtf_path,out_file):
    id2pos = defaultdict(list)
    with open(gtf_path,"r") as gtf:
        for line in gtf.readlines():
            if not line.startswith("#"):
                line = line.strip().split("\t")
                chr,type, start, end, score, strand, phase, attributes = line[0],line[2], line[3], line[4], line[5], line[6], line[
                    7], line[8]
                if type == "gene":
                    gene_id = re.search(r'gene_id \"(.+?)\";', attributes).group(1)
                    id2pos[gene_id] = [chr,start, end]
    print(len(id2pos))
    with open(out_file,"w") as w:
        for i in id2pos:
            w.write(i + "\t" + "\t".join(id2pos[i]) + "\n")

def extract_all_gene_pos(gtf_path,out_file):
    id2pos = defaultdict(list)
    id2locate = defaultdict(list)
    id2chr = defaultdict(str)
    with open(gtf_path,"r") as gtf:
        for line in gtf.readlines():
            if not line.startswith("#"):
                line = line.strip().split("\t")
                chr,type, start, end, score, strand, phase, attributes = line[0],line[2], line[3], line[4], line[5], line[6], line[
                    7], line[8]
                gene_id = re.search(r'gene_id \"(.+?)\";', attributes).group(1)
                id2pos[gene_id].extend([start, end])
                id2chr[gene_id] = chr
    for i in id2pos:
        id2locate[i] = [min(id2pos[i]), max(id2pos[i])]
    print(len(id2pos))
    with open(out_file,"w") as w:
        for i in id2locate:
            w.write(i + "\t" +id2chr[i]+"\t" +"\t".join(id2locate[i]) + "\n")

def extract_all_chr_length(chr_list,stat_file, out_file):
    with open(chr_list,"r") as l,open(stat_file,"r") as s,open(out_file,"w") as f:
        chrlist=[]
        chr2length={}
        for i in l.readlines():
            if i.strip().split("\t")[-1].lower() == "chromosome":
                chrlist.append(i.strip().split("\t")[0])
        for i in s.readlines():
            chr2length[i.strip().split("\t")[0]] = i.strip().split("\t")[1]
        for i in chrlist:
            f.write(i+"\t"+chr2length[i]+"\n")

def gtf_check_old(gtf_path,out_gtf_file):
    def is_in_interval(cds_dict,exon_dict):
        result = False
        start_loc = cds_dict[0]
        end_loc = cds_dict[1]
        for exon_region in exon_dict:
            exon_start = exon_region[0]
            exon_end = exon_region[1]
            if start_loc >= exon_start and end_loc <= exon_end:
                result = True
        return result

    all_exon_dict = defaultdict(list)
    #all_cds_dict = defaultdict(list)
    with open(gtf_path, "r") as gtf:
        for line in gtf.readlines():
            if not line.startswith("#"):
                line = line.strip().split("\t")
                chr, type, start, end, score, strand, phase, attributes = line[0], line[2], line[3], line[4], line[5],line[6], line[7], line[8]
                gene_id = re.search(r'gene_id \"(.+?)\";', attributes).group(1)
                if type == "exon":
                    all_exon_dict[gene_id].append([start,end])

    with open(gtf_path,"r") as gtf,open(out_gtf_file,"w") as out:
        for line in gtf.readlines():
            out.write(line)
            if not line.startswith("#"):
                line = line.strip().split("\t")
                chr, type, start, end, score, strand, phase, attributes = line[0], line[2], line[3], line[4], line[5], \
                                                                          line[6], line[7], line[8]
                gene_id = re.search(r'gene_id \"(.+?)\";', attributes).group(1)
                if type == "CDS":
                    cds_loc= [start,end]
                    is_exist = is_in_interval(cds_loc,all_exon_dict[gene_id])
                    if not is_exist:
                        line[2] = "exon"
                        out.write("\t".join(line)+"\n")

def gtf_check(gtf_path,out_gtf_file):
    def is_in_interval(cds_dict,exon_dict):
        result = False
        start_loc = cds_dict[0]
        end_loc = cds_dict[1]
        for exon_region in exon_dict:
            exon_start = exon_region[0]
            exon_end = exon_region[1]
            if start_loc >= exon_start and end_loc <= exon_end:
                result = True
        return result

    all_exon_dict = defaultdict(list)
    #all_cds_dict = defaultdict(list)
    with open(gtf_path, "r") as gtf:
        for line in gtf.readlines():
            if not line.startswith("#"):
                line = line.strip().split("\t")
                chr, type, start, end, score, strand, phase, attributes = line[0], line[2], line[3], line[4], line[5],line[6], line[7], line[8]
                gene_id = re.search(r'gene_id \"(.+?)\";', attributes).group(1)
                if type == "exon":
                    all_exon_dict[gene_id].append([start,end])

    with open(gtf_path,"r") as gtf,open(out_gtf_file,"w") as out:
        for line in gtf.readlines():
            if  line.startswith("#"):
                out.write(line)
            if not line.startswith("#"):
                line = line.strip().split("\t")
                chr, type, start, end, score, strand, phase, attributes = line[0], line[2], line[3], line[4], line[5], \
                                                                          line[6], line[7], line[8]
                # attributes  = re.sub(r'gene_name \"(.*?)([\t ])(.*?)\"',r'gene_name "\1_\3"',attributes)
                # raw_gene_name = re.search(r'gene_name \"(.*?)([ \']+)(.*?)\"', attributes)
                # if raw_gene_name:
                #     con = re.sub(r'[ \']', "_", raw_gene_name.group(2))
                #     new_name = "{}{}{}".format(raw_gene_name.group(1), con, raw_gene_name.group(3))
                #     attributes = re.sub(r'gene_name \"(.*?)([ \']+)(.*?)\"', r'gene_name "{}"'.format(new_name),
                #                         attributes)

                raw_gene_name = re.search(r'gene_name \"(.*?)\"', attributes)
                if raw_gene_name:
                    con = re.sub(r'[ \']', "_", raw_gene_name.group(1))
                    attributes = re.sub(r'gene_name \"(.*?)\"', r'gene_name "{}"'.format(con), attributes)
                line[8] = attributes
                out.write("\t".join(line) + "\n")
                gene_id = re.search(r'gene_id \"(.+?)\";', attributes).group(1)
                if type == "CDS":
                    cds_loc= [start,end]
                    is_exist = is_in_interval(cds_loc,all_exon_dict[gene_id])
                    if not is_exist:
                        line[2] = "exon"
                        out.write("\t".join(line)+"\n")



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='extract gene pos from gtf')
    # parser.add_argument('-g', '--genelist', dest='genelist', required=True, help='input genelist for pos extract')
    parser.add_argument('-r', '--refgtf', dest='refgtf', required=True, help='input raw gtf file path')
    parser.add_argument('-o', '--output', dest='output', required=True, help='output gene pos table')
    args = parser.parse_args()
    gtf_path = args.refgtf
    # gene_list = args.genelist
    out_file = args.output
    # extract_pos(gene_list,gtf_path,out_file)
    gtf_check(gtf_path,out_file)