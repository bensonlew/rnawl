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
import traceback

def extract_seqs(raw_data = None,ids = None ,level = None,type = None,output_dir = None):
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
        gene_list, trans_list, cds_pep_id_list, gene_dict = get_gene_details(raw_data, ids)
    elif level.lower() == "t":
        trans_list, cds_pep_id_list = get_trans_details(raw_data, ids)
    else:
        raise Exception("level must be in g or t")
    if level.lower() == "g":
        extract_gene_seqs(raw_data, gene_list, output_dir, gene_dict)
        if type:
            if "transcript" in type:
                extract_transcript_seqs(raw_data, trans_list, output_dir)
            if "cds" in type:
                extract_cds_seqs(raw_data, cds_pep_id_list, output_dir)
            if "pep" in type:
                extract_pep_seqs(raw_data, cds_pep_id_list, output_dir)
    elif level.lower() == "t":
        extract_transcript_seqs(raw_data, trans_list, output_dir)
        if type:
            if "cds" in type:
                extract_cds_seqs(raw_data, cds_pep_id_list, output_dir)
            if "pep" in type:
                extract_pep_seqs(raw_data, cds_pep_id_list, output_dir)

def extract_gene_seqs(raw_data, gene_list, output_dir, gene_dict):
        with open(os.path.join(raw_data, "txpt_seq"), "r") as f:
            gene_seq_dict = pickle.load(f)
        with open(os.path.join(output_dir, "extract_gene_seqs"), "w") as f:
            for gene_id in gene_list:
                try:
                    f.write(">{}".format(gene_dict[gene_id]) + "\n")
                    f.write(gene_seq_dict[gene_id] + "\n")
                except:
                    pass

def extract_transcript_seqs(raw_data, trans_list, output_dir):
        with open(os.path.join(raw_data, "txpt_seq"), "r") as f:
            trans_seq_dict = pickle.load(f)
        with open(os.path.join(output_dir, "extract_trans_seqs"), "w") as f:
            for trans_id in trans_list:
                try:
                    f.write(">{} | gene_id {gene_id}".format(trans_id,gene_id = t2g_dict[trans_id]) + "\n")
                    f.write(trans_seq_dict[trans_id] + "\n")
                except:
                    pass


def extract_cds_seqs(raw_data, cds_id_list, output_dir):
    with open(os.path.join(raw_data, "cds_seq"), "r") as f:
        cds_seq_dict = pickle.load(f)
    with open(os.path.join(output_dir, "extract_cds_seqs"), "w") as f:
        for cds_id in cds_id_list:
            try:
                f.write(">{} | trans_id {trans_id} |gene_id {gene_id}".format(cds_id,trans_id = c2t_dict[cds_id],gene_id = t2g_dict[c2t_dict[cds_id]]) + "\n"+cds_seq_dict[cds_id]+"\n")
            except Exception as e:
                print("error : {}".format(e))
                print("cds error :{}".format(cds_id))
                pass

def extract_pep_seqs(raw_data, pep_id_list, output_dir):
        with open(os.path.join(raw_data, "pep_seq"), "r") as f:
            pep_seq_dict = pickle.load(f)
        with open(os.path.join(output_dir, "extract_pep_seqs"), "w") as f:
            for pep_id in pep_id_list:
                try:
                    f.write(">{} cds_id {cds_id} |trans_id {trans_id} |gene_id {gene_id}".format(pep_id,
                                                                            cds_id = p2c_dict[pep_id],
                                                                            trans_id = c2t_dict[ p2c_dict[pep_id]],
                                                                            gene_id = t2g_dict[c2t_dict[ p2c_dict[pep_id]]]
                                                                            ) + "\n"+pep_seq_dict[pep_id] + "\n")
                except:
                    pass

def get_gene_details(raw_data,ids):
    trans_list = list()
    cds_pep_id_list = list()
    gene_list_new = list()
    gene_dict = dict()
    with open(os.path.join(raw_data,"gene_detail"),"r") as f:
        gene_detail_dict = pickle.load(f)
    with open(ids, "r") as g:
        all_gene = [line.strip() for line in g.readlines()]
    gene_list = set(all_gene)
    for gene_id in gene_list:
        gene_list_new.append(gene_detail_dict[gene_id]["gene_trans_id"])
        gene_dict[gene_detail_dict[gene_id]["gene_trans_id"]] = gene_id
        try:
            trans_ids = gene_detail_dict[gene_id]["tran_id"]
            trans_list.extend(trans_ids)
        except:
            pass
        try:
            cds_ids_list = gene_detail_dict[gene_id]["cds_pep_id"]
            for cds_id in cds_ids_list:
                if cds_id == '':
                    continue
                if ';' in cds_id:
                    cds_pep_id_list.extend(cds_id.strip().split(';'))
                else:
                    cds_pep_id_list.append(cds_id)
        except:
            pass

    return set(gene_list_new), set(trans_list), set(cds_pep_id_list), gene_dict

def get_trans_details(raw_data, ids):
    cds_pep_id_list = list()
    pep_id_list = list()
    with open(os.path.join(raw_data, "tran_detail"), "r") as f:
        trans_detail_dict = pickle.load(f)
    with open(ids, "r") as g:
        all_trans = [line.strip() for line in g.readlines()]
    trans_list = set(all_trans)
    for trans_id in trans_list:
        try:
            cds_ids = trans_detail_dict[trans_id]["cds_pep_id"]
            if cds_ids == '':
                continue
            if ';' in  cds_ids:
                cds_pep_id_list.extend(cds_ids.strip().split(';'))
            else:
                cds_pep_id_list.append(cds_ids)
        except:
            pass

    return set(trans_list), set(cds_pep_id_list)


def get_relation_dict(raw_data=None):
    with open(os.path.join(raw_data,"gene_detail"),"r") as f:
        t2g_dict = {}
        c2t_dict = {}
        p2c_dict = {}
        gene_detail_dict =pickle.load(f)
        for gene_id in gene_detail_dict:
            all_trans_ids = gene_detail_dict[gene_id]["tran_id"]
            for trans_id in all_trans_ids:
                t2g_dict[trans_id] = gene_id
    with open(os.path.join(raw_data, "tran_detail"), "r") as f:
        trans_detail_dict = pickle.load(f)
        for trans_id in trans_detail_dict:
            if trans_detail_dict[trans_id]["cds_pep_id"]:
                cds_ids = trans_detail_dict[trans_id]["cds_pep_id"]
                if ';' in cds_ids:
                    for cds_id in cds_ids.strip().split(';'):
                        c2t_dict[cds_id] = trans_id
                        p2c_dict[cds_id] = cds_id
                else:
                    cds_id = cds_ids
                    c2t_dict[cds_id] = trans_id
                    p2c_dict[cds_id] = cds_id


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
    t2g_dict, c2t_dict, p2c_dict = get_relation_dict(raw_data=args.r)
    extract_seqs(raw_data = args.r,ids = args.i ,level = args.l,type = args.t,output_dir = args.o)

