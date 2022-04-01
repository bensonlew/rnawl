# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import pandas as pd
import os
import re
import shutil
from Bio import SeqIO
from Bio import Phylo


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

def merge_16s(dir, out):
    files = os.listdir(dir)
    list1 = []
    with open(out, "w") as w:
        for file in files:
            name = file.split("_16S")[0]
            if os.path.getsize(dir + "/" + file) > 0:
                uniprot_iterator = list(SeqIO.parse(dir + "/" + file, "fasta"))
                if len(uniprot_iterator[0].seq) >= 800:
                    list1.append(uniprot_iterator[0].id)
                    w.write(">{}\n{}\n".format(name, uniprot_iterator[0].seq))
    return len(list1)



def add_merge(dir, type, out, sample_add=None, add_del=None):
    """
    将目录下的文件相同type的合并
    :param dir: 注释总目录
    :param type: 文件后缀类型，进行合并
    :param out: 输出文件
    :param sample_add: 是否在文件中加入样品名称
    :return:
    """
    files = os.listdir(dir)
    n = 1
    if os.path.exists(out):
        os.remove(out)
    for file in files:
        if re.search("{}".format(type), file):
            a = pd.read_table(dir + "/" + file, sep='\t', header=0, dtype={'fullVisitorId': 'str'})
            if sample_add in ['true', 'True', "TRUE"]:
                a['sample'] = file
            elif sample_add in ['false', 'False', "FALSE"]:
                pass
            if add_del:
                a = a.drop([add_del], axis=1)
            if n == 1:
                a.to_csv(out, mode='a', sep='\t', header=True, index=False)
            elif n > 1:
                a.to_csv(out, mode='a', sep='\t', header=0, index=False)
            n += 1


def add_fasta(dir,type, out):
    """
    将目录下的文件进行合并，
    :param type:文件后缀类型+-
    :param dir:处理文件目录
    :param out:处理完后结果目录
    :return:
    """
    files = os.listdir(dir)
    if os.path.exists(out):
        os.remove(out)
    list = []
    for file in files:
        if re.search("{}".format(type), file):
            list.append(dir + "/" + file)
    des = " ".join(list)
    os.system("cat {} > {}".format(des, out))

def get_fasta(file, sampl_list, out):
    uniprot_iterator = SeqIO.parse(file, "fasta")
    records = list(uniprot_iterator)
    with open (out, "w") as g:
        for sample in sampl_list:
            for t in records:
                if sample in t.id:
                    g.write(">{}\n{}\n".format(sample, t.seq))


def vfdb_anno(core_file, predict_file, all_file):
    with open (core_file, "r") as f1,open (predict_file, "r") as f2,open (all_file, "w") as g:
        lines = f1.readlines()
        names = lines[0].strip().split("\t")
        names.append("Type")
        g.write("\t".join(names)+ "\n")
        for line in lines[1:]:
            lin = line.strip().split("\t")
            lin.append("Core")
            g.write("\t".join(lin)+ "\n")
        lines2 = f2.readlines()
        for line in lines2[1:]:
            lin = line.strip().split("\t")
            lin.append("Predict")
            g.write("\t".join(lin) + "\n")


def get_sample_from_tree(file):
    list =[]
    tree = Phylo.read(file, "newick")
    for i in tree.get_terminals():
        list.append(str(i))
    return ",".join(list)


def get_gff_files(dir, outdir):
    for sample in os.listdir(dir):
        if os.path.isdir(dir + "/" + sample):
            with open(outdir + "/" + sample + ".gff", "w") as g:
                g.write("Gene ID\tSequence id\tStart\tEnd\tStrand\tGene Length(bp)\tProtein Length\tA.start\tA.end\tInitiator Codon\tTerminator Codon\tType\n")
                for file in os.listdir(dir + "/" + sample):
                    if file.endswith(("_CDS.gff")) and get_num(dir + "/" + sample + "/" + sample + "_CDS.gff") > 1:
                        with open(dir + "/" + sample + "/" + sample + "_CDS.gff", "r") as f:
                            lines = f.readlines()
                            for line in lines[1:]:
                                lin = line.strip().split("\t")
                                g.write("{}\tgene\n".format("\t".join(lin)))
                    elif file.endswith(("_tRNA.gff")) and get_num(dir + "/" + sample + "/" + sample + "_tRNA.gff") > 1:
                        with open(dir + "/" + sample + "/" + sample + "_tRNA.gff", "r") as f:
                            lines = f.readlines()
                            for line in lines[1:]:
                                lin = line.strip().split("\t")
                                location = lin[1].split("_tRNA")[0] + "_ORF_tRNA"
                                if int(lin[3]) > int(lin[2]):
                                    strand = "+"
                                    len = int(lin[3]) - int(lin[2]) + 1
                                elif int(lin[2]) > int(lin[3]):
                                    len = int(lin[2]) - int(lin[3]) + 1
                                    strand = "-"
                                des = [lin[0], location, str(lin[2]), str(lin[3]), strand, str(len), "-", "-", "-", "-",
                                       "-"]
                                g.write("{}\ttRNA\n".format("\t".join(des)))
                    elif file.endswith(("_rRNA.gff")) and get_num(dir + "/" + sample + "/" + sample + "_rRNA.gff") > 1:
                        with open(dir + "/" + sample + "/" + sample + "_rRNA.gff", "r") as f:
                            lines = f.readlines()
                            for line in lines[1:]:
                                lin = line.strip().split("\t")
                                location = lin[1].split("_rRNA")[0] + "_ORF_rRNA"
                                if lin[5] in ["+"]:
                                    len = int(lin[3]) - int(lin[2]) + 1
                                elif lin[5] in ["-"]:
                                    len = int(lin[2]) - int(lin[3]) + 1
                                des = [lin[0], location, lin[2], lin[3], lin[5], str(len), "-", lin[8], lin[9], "-",
                                       "-"]
                                g.write("{}\trRNA\n".format("\t".join(des)))


def get_num(file):
    with open(file, "r") as f:
        lines = f.readline()
    return len(lines)

def sort_seq(input, out_seq):
    dict = {}
    for i in SeqIO.parse(input, "fasta"):
        dict[i] = len(i.seq)
    dict2 = sorted(dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    list = []
    for i in dict2:
        list.append(i[0])
    SeqIO.write(list, out_seq, "fasta")


def format_category(is_percent,origin_category,max_sample_number,category_list):
    """
    对传入的category进行format，有两种形式：数值型和百分比
    数值型：{'core'：{"min":"N", "max": "N"}, 'dispensable'：{"min":"2", "max": "N-1"}, 'unique'{"min":"0", "max": "1"}}
    百分比：{core : {min: 0.95, max: 1}, dis:{min :0.05 , max : 0.95}, unique : {min:0, max:0.05}}
    :param is_percent: 是否是百分比
    :param origin_category: 原始的名称
    :param max_sample_number: 当前分析的样本数
    :param category_list: 分类方案内容
    :return:
    add by qingchen.zhang @20200917
    """
    category_dict = {}
    if is_percent: #根据传入的是百分比还是整数进行判断
        for categ in category_list:
            categ_dict = {}
            categ_dict["min"] = float(max_sample_number * float(origin_category[categ]["min"]))
            categ_dict["max"] = float(max_sample_number * float(origin_category[categ]["max"]))
            category_dict[categ] = categ_dict
    else:
        for categ in category_list:
            categ_dict = {}
            if categ in origin_category.keys():## 增加判断，根据后面的需求将样本数替换为N
                if re.search(r"N", str(origin_category[categ]["min"])):
                    if re.search(r"-", str(origin_category[categ]["min"])):
                        min_nu = max_sample_number - int(str(origin_category[categ]["min"]).split("-")[1])
                    else:
                        min_nu = max_sample_number
                else:
                    min_nu = int(float(origin_category[categ]["min"]))

                if re.search(r"N", str(origin_category[categ]["max"])):
                    if re.search(r"-", str(origin_category[categ]["max"])):
                        max_nu = max_sample_number - int(str(origin_category[categ]["max"]).split("-")[1])
                    else:
                        max_nu = max_sample_number
                else:
                    max_nu = int(float(origin_category[categ]["max"]))
                categ_dict["min"] = int(float(min_nu))
                categ_dict["max"] = int(float(max_nu))
            category_dict[categ] = categ_dict
    return category_dict


def add_signature(cluster_data, category,category_list, is_percent, max_sample_number):
    """
    目的：为cluster_data中的每一条cluster加标签
    :param cluster_data: dataframe
    :param category: 分类标签
    :param category_list: 分类方案
    :param is_percent: 是否是百分比
    :param max_sample_number: 当前分析的样本数
    add by qingchen.zhang @20200917
    :return:
    """
    for x in cluster_data.index: #对每个cluster加标签如：core、dispensable、unique
        for categ in category_list:
            sample_number = int(cluster_data.loc[x, "Sample_number"])
            gene_number = int(cluster_data.loc[x, "Gene_number"])
            if is_percent: #根据传入的是百分比进行判断
                if category[categ]["min"] == category[categ]["max"]:
                    if (sample_number == max_sample_number) and (category[categ]["min"] >= sample_number):##core样本
                        cluster_data.loc[x, "Category"] = categ
                    elif (sample_number == 1) and (category[categ]["min"] <= sample_number):##unique样本
                        cluster_data.loc[x, "Category"] = categ
                else:
                    if (category[categ]["min"] <= sample_number) and (sample_number < category[categ]["max"]): ## dispensiable 这里取小不取大
                        cluster_data.loc[x, "Category"] = categ
                    elif (sample_number == max_sample_number) and (category[categ]["max"] == sample_number) and (category[categ]["min"] <= sample_number) and (categ in ["core", "Core", "CORE"]):##core样本，如果选择的core基因的名称不为core可能会存在问题
                        cluster_data.loc[x, "Category"] = categ
            else:#根据是数值进行判断
                if category[categ]["min"] == category[categ]["max"]:
                    if (sample_number == max_sample_number) and (category[categ]["min"] >= sample_number):##core样本
                        cluster_data.loc[x, "Category"] = categ
                    elif (sample_number == 1) and (category[categ]["min"] <= sample_number):##unique样本
                        cluster_data.loc[x, "Category"] = categ
                else:
                    if (category[categ]["min"] <= sample_number) and (sample_number <= category[categ]["max"]): ## dispensiable 这里都取
                        cluster_data.loc[x, "Category"] = categ
                    elif (sample_number == max_sample_number) and (category[categ]["max"] == sample_number) and (category[categ]["min"] <= sample_number) and (categ in ["core", "Core", "CORE"]):##core样本
                        cluster_data.loc[x, "Category"] = categ
    return cluster_data

def chang_value(input, output, type=None):
    data = pd.read_table(input, sep='\t', header=0)
    a = list(data.columns)
    del (a[0])
    if type:
        type = type
    else:
        type = 100
    for i in a:
        data[i] = data[i].apply(lambda x: type - x)
    data.to_csv(output, sep='\t', header=True, index=False)

def anno_kegg(list, output):
    for i in list:
        sample = os.path.basename(i).split(".kegg_anno")[0]
        if list.index(i) == 0:
            data = pd.read_table(i, sep='\t', header=0)
            data.rename(columns={'#Query': 'Gene ID'}, inplace=True)
            data['sample'] = sample
            data.to_csv(output, mode='a', sep='\t', header=True, index=False)
        else:
            data = pd.read_table(i, sep='\t', header=0)
            data.rename(columns={'#Query': 'Gene ID'}, inplace=True)
            data['sample'] = sample
            data.to_csv(output, mode='a', sep='\t', header=0, index=False)
