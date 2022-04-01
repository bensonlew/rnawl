#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2020/5/11 16:31
@file    : meta_Norm.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""

import multiprocessing
import pandas as pd
import numpy as np
import re
import argparse
from mako.template import Template
import os





def run_norm(args, type_='pos'):
    dir_ = './'  # os.path.dirname(args.pos)

    dir_ = os.path.join(dir_, type_)
    if not os.path.exists(dir_):
         os.mkdir(dir_)
    os.chdir(dir_)
    ord_file = args.order
    smp_list = []
    smp_list_new = []

    with open(ord_file, 'r') as ord_r, open("sample.info.csv",'w') as smp_w:
        ord_r.readline()
        smp_w.write("sample.name,injection.order,class\n")
        for line in ord_r:
            items = line.strip().split('\t')
            smp_list.append(items[0])
            smp_name = re.sub(r'[-_]', '.', items[0])
            smp_list_new.append(smp_name)
            if 'QC' in smp_name:
                mod = "QC"
            else:
                mod = "Subject"
            smp_w.write("%s,%s,%s\n" % (smp_name, items[1], mod))

    sample_dict = dict(zip(smp_list_new, smp_list))
    if type_ == 'pos':
        in_file = args.pos
    else:
        in_file = args.neg
    try:
        with open(in_file, 'rb') as ir:
            dat_df = pd.read_excel(ir, sep='\t')
    except:
        with open(in_file, 'rb') as ir:
            dat_df = pd.read_csv(ir, sep='\t', quoting=0, quotechar='"', encoding='utf-8')
    del_id = False
    if args.mark not in dat_df.columns:
        del_id = True
        dat_df.insert(0, args.mark, np.arange(0, dat_df.shape[0]) + 1)
    header = dat_df.columns.tolist()
    peak_name = args.mark
    mayheads = ['m/z', 'M/Z', 'rt', 'RT (min)', 'Retention time', 'mz']

    chs_head = [x for x in mayheads if x in header]
    chs_head.insert(0, peak_name)
    new_header = ["name", "mz", "rt"]
    chs_head.extend(smp_list)
    deal_dat = dat_df[chs_head]
    new_header.extend(smp_list_new)
    renm_head = dict(zip(chs_head, new_header))
    deal_dat.rename(columns=renm_head, inplace=True)
    with open("data.csv", 'w') as dw:
        deal_dat.to_csv(dw, sep=",", index=False)
    des_header = [x for x in header if x not in smp_list]
    des_df = dat_df[des_header]

    Rcmd=r"""
library(dplyr)
library(ggplot2)
library(patchwork)
source("${package_path}/metabolome/MetNormalizer/checkData.R")
source("${package_path}/metabolome/MetNormalizer/metNor.R")
source("${package_path}/metabolome/MetNormalizer/SXTsvrNor.R")
source("${package_path}/metabolome/MetNormalizer/tools.R")
#metNor(ms1.data.name = "data.csv",sample.info.name = "sample.info.csv",minfrac.qc = as.numeric(qc_f),minfrac.sample = as.numeric(smp_f), multiple = as.numeric(top))
# library(MetNormalizer)
# library(snow)
print(1111111)
qc_f <- '${qc_f}'
smp_f <- '${smp_f}'
method <- '${method}'
top <- '${top}'
#MetNormalizer(ms1.data.name = "data.csv",sample.info.name = "sample.info.csv",minfrac.qc = as.numeric(qc_f),minfrac.sample = as.numeric(smp_f), normalization.method = method,multiple = as.numeric(top))
metNor(ms1.data.name = "data.csv",sample.info.name = "sample.info.csv",minfrac.qc = as.numeric(qc_f),minfrac.sample = as.numeric(smp_f), multiple = as.numeric(top))
    """
    f = Template(Rcmd)
    Rcmd_ = f.render(qc_f=args.qc_f, smp_f=args.smp_f, method=args.method, top=args.top, package_path=args.rsrc_path)
    #global r_cmd
    #Rcmd_ = r_cmd + '\n' + Rcmd_
    with open('meta_norm.r', 'w') as mw:
        mw.write(Rcmd_)
    print(Rcmd_)
    # robj.r(Rcmd_)
    # robj.r.source('meta_norm.r')
    # os.system('start "meta norm" /MAX D:\\R_box\\R-3.6.1\\bin\\Rscript.exe meta_norm.r')
    os.system('Rscript  meta_norm.r')

    # Dealed_data = pd.read_csv('svr normalization result\\data svr nor.csv', sep=",", index_col=0)
    Dealed_data = pd.read_csv('svr_normalization_result/data_svr_normalization.csv', sep=",", index_col=0)
    drop_col = ["mz", "rt", "sample.nor.rsd", "QC.nor.rsd"]
    Dld_hds = Dealed_data.columns.tolist()
    new_hds = [x for x in Dld_hds if x not in drop_col]
    Dld_data_Dld = Dealed_data[new_hds]
    Dld_data_Dld = Dld_data_Dld[Dld_data_Dld.sum(axis=1) > 0]

    result = pd.merge(des_df, Dld_data_Dld, how="inner", left_on=peak_name, right_on='name').drop('name', axis=1)
    result.rename(columns=sample_dict, inplace=True)
    if del_id:
        result = result.drop(args.mark, axis=1)
    result.to_csv("Dealed_%s" % os.path.basename(in_file), sep="\t", index=False, quotechar='"', quoting=3, encoding='utf-8', header=True)
    os.chdir('../')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Normalization and Integration of Large-Scale Metabolomics Data Using Support Vector Regression")
    parser.add_argument("-pos", type=str, help="data matrix needed to be dealed")
    parser.add_argument("-neg", type=str, help="data matrix needed to be dealed")
    parser.add_argument("-mark", type=str,  help="the name of the peak name column,'ID' or 'Metabolite'")
    parser.add_argument('-order', type=str,  help="the file record injection order of samples")
    parser.add_argument("-qc_f", type=str,
                        help="If the zero value rate of one peaks in all the QC samples is bigger than 1- minfrac.qc, this peak will be removed from the dataset. Default is 0, this means no peaks will be removed")
    parser.add_argument("-smp_f", type=str,
                        help="If the zero value rate of one peaks in all the subject samples is bigger than 1- minfrac.sample, this peak will be removed from the dataset. Default is 0, this means no peaks will be removed")
    parser.add_argument("-method", type=str,  help="normalization.method 'svr','loess' or 'all'")
    parser.add_argument("-top", type=str,
                        help="using top mutiple peaks correlated peaks,for'tof'data,top=5,else top =1")
    parser.add_argument("-rsrc_path", type=str,
                        help="metNormalizer path")


    args = parser.parse_args()


    if args.neg:
        # process_p = multiprocessing.Process(target=run_norm, args=(args, 'pos',))
        # process_n = multiprocessing.Process(target=run_norm, args=(args, 'neg',))
        # process_p.start()
        # process_n.start()
        # process_p.join()
        # process_n.join()
        run_norm(args, 'pos')
        run_norm(args, 'neg')

    else:
        run_norm(args, 'pos')
