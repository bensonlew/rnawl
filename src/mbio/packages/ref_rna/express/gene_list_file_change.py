#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import re
from itertools import islice 
"""对表达量和差异分析的结果进行修改，符合输出给客户的标准"""

def check_diff_fpkm(input_diff, out_diff_name, fpkm, count):
    """
    :params input_diff: 'gene:ENSEG' 格式的差异分析输出文件
    :params out_diff_name: 输出符合差异分析的文件名
    :params fpkm: 基因的fpkm表
    :params count: 基因的count表
    """
    input_params=[input_diff, fpkm, count]
    check_file=map(os.path.exists, input_params)
    print check_file
    #for i,index in enumerate(check_file):
    #    if not i:
    #        raise Exception as e:
    #            print "{}".format(input_params[index])+"不存在！"
    out_diff=os.path.join(os.path.split(input_diff)[0],out_diff_name)
    with open(input_diff,"r+") as f1, open(out_diff,"w+") as f2, open(fpkm,"r+") as f3, open(count,"r+") as f4:
        line1=islice(f1,1,None)
        line3=islice(f3,1,None)
        line4=islice(f4,1,None)
        data=[i.strip().split("\t")for i in line1]

        _fpkm=[i.strip().split("\t")for i in line3]
        fpkm_data={}
        for i in _fpkm:
            fpkm_data[i[0]]=i
        _count=[i.strip().split("\t")for i in line4]
        count_data={}
        for i in _count:
            count_data[i[0]]=i
        
        for j in data:
            if re.search("gene:", j[0]):
                geneid=re.sub('gene:',"",j[0])
                j.pop(0)
                j.insert(0,geneid)
            else:
                geneid=j[0]
                
            _line=geneid
            if count_data.has_keys(geneid):
                count_value=count_data[geneid][1:]
                _line+="\t".join(count_value)
                if fpkm_data.has_keys(geneid):
                    fpkm_value=fpkm_data[geneid][1:]
                    _line+="\t".join(fpkm_value)+"\t".join(j[1:])+"\n"
                    f2.write(_line)
    return out_diff
"""输出全部基因列表"""
def exp_gene_list(count_path):
    """
    :params count_path: 基因的count表
    """
    if os.path.exists(count_path):
        pass
    gene_list=os.path.join(os.path.split(count_path)[0],"all_gene_list.txt")
    with open(count_path,'r+') as f,open(gene_list,'w+') as w:
        line1=islice(f,1,None)
        for line in line1:
            line.strip()
            w.write(line.split("\t")[0]+"\n")
    return gene_list
    
