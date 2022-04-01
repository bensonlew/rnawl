# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
#last_modified:20161101 by qindanhua

import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import kendalltau
import numpy as np
import os

#input files: otu table/species table and factor table locations
def pearsonsCorrelation(out_species_location,factor_location,exportfilepath1,exportfilepath2,method="pearsonr"):
    if not os.path.isfile(out_species_location):
        raise Exception('OTU或Species文件不存在:{}'.format(out_species_location))
    if not os.path.isfile(factor_location):
        raise Exception('样本环境文件不存在:{}'.format(factor_location))
    pvaluefile=open(exportfilepath2,'w')
    with open(exportfilepath1,"w") as newfile:
        df1=pd.read_csv(out_species_location,sep="\t", dtype=str)
        df2=pd.read_csv(factor_location,sep="\t", dtype=str)
        #get otu ids by column name
        otuids=list(df1[df1.columns[0]])
        '''
        try:
            otuids=list(df1['OTUID'])
        except KeyError:
            otuids=list(df1['OTU'])
        '''
        #get sample ids from otu table
        analyzed_ids=list(df1.columns[1:])
        #get whole sample ids from factor table
        all_factors=list(df2.columns)[1:]
        column1=list(df2.columns)[0]
        allsamples={}
        k=0
        for samplename in list(df2[column1]):
            allsamples[samplename]=k
            k+=1
        newfile.write('name'+'\t')
        pvaluefile.write('name'+'\t')
        for j in range(0,len(all_factors)):
            pvaluefile.write(all_factors[j])
            newfile.write(all_factors[j])
            if j!=len(all_factors)-1:
                newfile.write('\t')
                pvaluefile.write('\t')
        newfile.write('\n')
        pvaluefile.write('\n')
        for i in range(0,len(otuids)):
            otu = df1.loc[i][0].split("; ")[-1]  # modify by zhujuan for function name
            #otu=df1.loc[i][0]  # modify by zhujuan for function name
            corr_write_line = [otu]
            pavalue_write_line = [otu]
            isnan = False
            first_array=np.asarray(map(lambda x:float(x),df1.loc[i][1:]))
            # print first_array
            for f in all_factors:
                degree=[]
                for item in analyzed_ids:
                    degree.append(float(df2.loc[allsamples[item],f]))
                second_array=np.asarray(degree)
                # print second_array
                if method == "pearsonr":
                    correlation=pearsonr(first_array,second_array)
                if method == "spearmanr":
                    correlation=spearmanr(first_array,second_array)
                if method == "kendalltau":    ## add for kendall by zhujuan 2017.08.30
                    correlation=kendalltau(first_array,second_array)
                if np.isnan(correlation[0]):
                    isnan = True
                corr_write_line.append(str(correlation[0]))
                pavalue_write_line.append(str(correlation[1]))
            if isnan:
                continue
            else:
                newfile.write("\t".join(corr_write_line)+"\n")
                pvaluefile.write("\t".join(pavalue_write_line)+"\n")

import sys
pearsonsCorrelation(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4], sys.argv[5])
