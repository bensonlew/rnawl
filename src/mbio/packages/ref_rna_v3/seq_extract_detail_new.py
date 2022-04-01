# coding=utf-8

from Bio import SeqIO
import fileinput
import re
import os
import subprocess
import urllib2
from collections import defaultdict
import regex
# import pandas as pd
# from collections import Counter
# import matplotlib.pyplot as plt
import sys
import pickle

def extract_seqs(raw_data = None,ids = None ,level = None,type = None,output_dir = None):
# def extract_seqs(raw_data = None,ids = None ,level = None,type = None,output_dir = None,t2g= None,c2t = None,p2c = None):
    """
    步长统计
    :param fasta_file: 输入的fa文件
    :param fasta_to_txt:输出的统计数据txt
    :param group_num:按照步长统计几组
    :param step:步长
    :param stat_out:统计的数据汇总信息txt（fasta_to_txt文件的汇总）
    :return:
    """
    if level.lower() == "g":
        gene_list, trans_list, cds_id_list, pep_id_list = get_gene_details(raw_data,ids)
    elif level.lower() == "t":
        trans_list, cds_id_list, pep_id_list = get_trans_details(raw_data, ids)
    else:
        raise Exception("level must be in g or t")
    if level.lower() == "g":
        extract_gene_seqs(raw_data,gene_list,output_dir)
        if type:
            if "transcript" in type:
                extract_transcript_seqs(raw_data,trans_list,output_dir)
            if "cds" in type:
                extract_cds_seqs(raw_data,cds_id_list,output_dir)
            if "pep" in type:
                extract_pep_seqs(raw_data,pep_id_list,output_dir)
    elif level.lower() == "t":
        extract_transcript_seqs(raw_data, trans_list, output_dir)
        if type:
            if "cds" in type:
                extract_cds_seqs(raw_data,cds_id_list,output_dir)
            if "pep" in type:
                extract_pep_seqs(raw_data,pep_id_list,output_dir)

def extract_gene_seqs(raw_data,gene_list,output_dir):
        with open(os.path.join(raw_data,"gene_seq"),"r") as f:
            gene_seq_dict = pickle.load(f)
        with open(os.path.join(output_dir,"extract_gene_seqs"),"w") as f:
            for gene_id in gene_list:
                try:
                    f.write(">{}".format(gene_id)+"\n")
                    f.write(gene_seq_dict[gene_id]+"\n")
                except:
                    pass

def extract_transcript_seqs(raw_data,trans_list,output_dir):
        with open(os.path.join(raw_data,"txpt_seq"),"r") as f:
            trans_seq_dict = pickle.load(f)
        with open(os.path.join(output_dir,"extract_trans_seqs"),"w") as f:
            for trans_id in trans_list:
                try:
                    if t2g_dict:
                        f.write(">{} | gene_id {}".format(trans_id,t2g_dict[trans_id])+"\n")
                    else:
                        f.write(">{}".format(trans_id) + "\n")
                    f.write(trans_seq_dict[trans_id]+"\n")
                except:
                    pass


def extract_cds_seqs(raw_data, cds_id_list, output_dir):
    with open(os.path.join(raw_data, "cds_seq"), "r") as f:
        cds_seq_dict = pickle.load(f)
    with open(os.path.join(output_dir, "extract_cds_seqs"), "w") as f:
        for cds_id in cds_id_list:
            if c2t_dict:
                try:
                    f.write(">{} |trans_id {} |gene_id {}".format(cds_id,c2t_dict[cds_id],t2g_dict[c2t_dict[cds_id]]) + "\n"+cds_seq_dict[cds_id]+"\n")
                except:
                    try:
                        f.write(">{} |trans_id {} |gene_id {}".format(cds_id,c2t_dict[cds_id],t2g_dict[c2t_dict[cds_id]]) + "\n"+cds_seq_dict[cds_id]["sequence"] + "\n")
                    except:
                        pass
            else:
                try:
                    f.write(">{}".format(cds_id) + "\n"+cds_seq_dict[cds_id]+"\n")
                except:
                    try:
                        f.write(">{}".format(cds_id) + "\n"+cds_seq_dict[cds_id]["sequence"] + "\n")
                    except:
                        pass

def extract_pep_seqs(raw_data, pep_id_list, output_dir):
        with open(os.path.join(raw_data, "pep_seq"), "r") as f:
            pep_seq_dict = pickle.load(f)
        with open(os.path.join(output_dir, "extract_pep_seqs"), "w") as f:
            for pep_id in pep_id_list:
                if p2c_dict:
                    try:
                        f.write(">{} |cds_id {} |trans_id {} |gene_id {}".format(pep_id,p2c_dict[pep_id],c2t_dict[p2c_dict[pep_id]],t2g_dict[c2t_dict[p2c_dict[pep_id]]]) + "\n"+pep_seq_dict[pep_id] + "\n")
                    except:
                        try:
                            f.write(">{} |cds_id {} |trans_id {} |gene_id {}".format(pep_id,p2c_dict[pep_id],c2t_dict[p2c_dict[pep_id]],t2g_dict[c2t_dict[p2c_dict[pep_id]]]) +"\n"+pep_seq_dict[pep_id]["sequence"] + "\n")
                        except:
                            pass
                else:
                    try:
                        f.write(">{}".format(pep_id) + "\n"+pep_seq_dict[pep_id] + "\n")
                    except:
                        try:
                            f.write(">{}".format(pep_id) + "\n"+pep_seq_dict[pep_id]["sequence"] + "\n")
                        except:
                            pass

def get_gene_details(raw_data,ids):
    gene_list =set()
    trans_list =set()
    cds_id_list =set()
    pep_id_list =set()
    with open(os.path.join(raw_data,"gene_detail"),"r") as f:
        gene_detail_dict =pickle.load(f)
    with open(ids,"r") as g:
        all_gene = [line.strip() for line in g.readlines()]
    gene_list =set(all_gene)
    for gene_id in gene_list:
        try:
            trans_ids = gene_detail_dict[gene_id]["trans_id"]
            trans_list = trans_list | set(trans_ids)
        except:
            pass
        try:
            cds_ids = gene_detail_dict[gene_id]["cds_id"]
            for cds_id in cds_ids:
                cds_id_list = cds_id_list | set(cds_id)
        except:
            pass
        try:
            pep_ids = gene_detail_dict[gene_id]["pep_id"]
            for pep_id in pep_ids:
                pep_id_list = pep_id_list | set(pep_id)
        except:
            pass
    return gene_list,trans_list,cds_id_list,pep_id_list

def get_trans_details(raw_data,ids):
    trans_list =set()
    cds_id_list =set()
    pep_id_list =set()
    with open(os.path.join(raw_data,"trans_detail"),"r") as f:
        trans_detail_dict =pickle.load(f)
    with open(ids,"r") as g:
        all_trans = [line.strip() for line in g.readlines()]
    trans_list =set(all_trans)
    for trans_id in trans_list:
        try:
            cds_ids = trans_detail_dict[trans_id]["cds_id"]
            cds_id_list = cds_id_list | set(cds_ids)
        except:
            pass
        try:
            pep_ids = trans_detail_dict[trans_id]["pep_id"]
            pep_id_list = pep_id_list | set(pep_ids)
        except:
            pass
    return trans_list,cds_id_list,pep_id_list


def get_relation_dict(raw_data=None):
    with open(os.path.join(raw_data,"gene_detail"),"r") as f:
        t2g_dict = {}
        c2t_dict = {}
        p2c_dict = {}
        gene_detail_dict =pickle.load(f)
        for gene_id in gene_detail_dict:
            all_trans_ids = gene_detail_dict[gene_id]["trans_id"]
            for trans_id in all_trans_ids:
                t2g_dict[trans_id] = gene_id
    with open(os.path.join(raw_data, "trans_detail"), "r") as f:
        trans_detail_dict = pickle.load(f)
        for trans_id in trans_detail_dict:
            if trans_detail_dict[trans_id]["cds_id"]:
                for n,cds_id in enumerate(trans_detail_dict[trans_id]["cds_id"]):
                        c2t_dict[cds_id] =  trans_id
                        pep_id = trans_detail_dict[trans_id]["pep_id"][n]
                        p2c_dict[pep_id] = cds_id
    return t2g_dict,c2t_dict,p2c_dict



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-r', type=str, metavar="raw_data", required=True,
                        help="please input fasta file ")
    parser.add_argument('-i', type=str, metavar="gene_transcript_ids", help="id list for extract ", required=True)
    parser.add_argument('-l', type=str, metavar="level", help="gene or transcript level", required=True)
    parser.add_argument('-o', type=str, metavar="output_dir",default=None, help="default is local dir. Output directory name", required=True)
    parser.add_argument('-t', type=str, metavar="extract_type", default=None, help="extract_types ")
    #
    args = parser.parse_args()
    t2g_dict,c2t_dict,p2c_dict = get_relation_dict(raw_data =  args.r)
    extract_seqs(raw_data=args.r, ids=args.i, level=args.l, type=args.t, output_dir=args.o)
    # extract_seqs(raw_data = args.r,ids = args.i ,level = args.l,type = args.t,output_dir = args.o ,t2g= t2g_dict,c2t = c2t_dict,p2c = p2c_dict)

