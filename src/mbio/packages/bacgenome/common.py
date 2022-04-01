# -*- coding: utf-8 -*-
# __author__ = 'ysh'
# last_modify:20190401
import pandas as pd
from Bio import SeqIO
import os
import shutil

def sum_stat(file_path,cname1,cname2,out,stat_method="sum",ocname=None):
    '''
    根据基因数统计stat表
    :param file_path: input file name
    :param cname1: columns name1
    :param cname2: columns name2
    :param out: outputfile name
    :param stat_method: use sum or count method
    :param ocname: outputfile second columns name if not given, name = "Num"
    '''
    table = pd.read_table(file_path, sep="\t",header=0)
    select_table = table[[cname1,cname2]]
    if stat_method == "sum":
        result = select_table.groupby(cname1).sum()
    elif stat_method == "count":
        result = select_table.groupby(cname1).count()
    if ocname:
        result.columns = [ocname]
    else:
        result.columns = ["Num"]
    result.to_csv(out,sep="\t")


def add_coverge(gene_anno, align_file):
    anno_table = pd.read_table(gene_anno, sep="\t",header=0)
    align_table = pd.read_table(align_file, sep="\t",header=0)


def fasta_cutoff(fa, cut_off, out):
    with open(out, "w") as w:
        uniprot_iterator = SeqIO.parse(fa, "fasta")
        records = list(uniprot_iterator)
        for i in records:
            if len(list(i.seq)) >= int(cut_off):
                w.write(">{}\n{}\n".format(i.id, i.seq))

def get_scaffoldlist(fa):
    list1 =[]
    for records in SeqIO.parse(fa, "fasta"):
        list1.append(records.id)
    return list1

def get_namelist(file):
    list1 =[]
    with open (file, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            name =line.strip().split("\t")[0]
            list1.append(name)
    return list(set(list1))

def get_num(file):
    with open(file, "r") as f:
        lines =f.readlines()
        num = len(lines)
    return num


def link_dir(olddir, newdir):
    """
    hard link directory from olddir to newdir
    :param olddir:
    :param newdir:
    :return:
    """
    if not os.path.isdir(olddir):
        raise Exception("不存在路径: %s" % olddir)
    allfiles = os.listdir(olddir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    newfiles = [os.path.join(newdir,i) for i in allfiles]
    for newfile in newfiles:
        if os.path.exists(newfile):
            if os.path.isfile(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile):
                shutil.rmtree(newfile)
    if len(allfiles) >= 1:
        for i in allfiles:
            if os.path.isfile(os.path.join(olddir, i)):
                os.link(os.path.join(olddir, i),os.path.join(newdir, i))
            elif os.path.isdir(os.path.join(olddir, i)):
                link_dir(os.path.join(olddir, i),os.path.join(newdir, i))
    else:
        raise Exception("结果文件夹为空:%s" % allfiles)

def link_file(oldfile, newfile):
    """
    hard link file from oldfile to newfile
    :param oldfile:
    :param newfile:
    :return:
    """
    if not os.path.isfile(oldfile):
        raise Exception("不存在文件：%s" % oldfile)
    if os.path.exists(newfile):
        os.remove(newfile)
    os.link(oldfile, newfile)




