# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import pandas as pd
import os
import re
import shutil
from Bio import SeqIO
from Bio import Entrez
from biocluster.config import Config


def get_num(file):
    with open(file, "r") as f:
        lines =f.readline()
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