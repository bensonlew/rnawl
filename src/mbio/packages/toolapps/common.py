# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import pandas as pd
import os
import re
import shutil
from Bio import SeqIO
from Bio import Entrez
from biocluster.config import Config
from collections import defaultdict


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


def get_info(dir):
    sample_path = defaultdict(list)
    with open(dir + "/list.txt") as fr:
        for line in fr:
            tmp = line.strip().split('\t')
            if tmp[1] in sample_path.keys():
                if tmp[2] == 'l':
                    sample_path[tmp[1]].insert(0, dir + '/' + tmp[0])
                else:
                    sample_path[tmp[1]].append(dir + '/' + tmp[0])
            else:
                sample_path[tmp[1]].append(dir + '/' + tmp[0])
    return sample_path

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