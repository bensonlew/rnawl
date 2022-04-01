#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = konghualei
# last_modify: 2016.11.15
import os
import shutil
"""根据样本的分组信息，转换样本名为sample1_1等格式，用于差异表达分析"""

def changesamplename(file_path_count, file_path_fpkm, out_count, out_fpkm, group_file, venn_table_path=None):
    """
    :params file_path_count, 输入表达量count文件
    :params file_path_fpkm, 输入表达量fpkm文件
    :params out_count, 生成的count表达量的文件名 例如"count.txt"
    :params out_fpkm, 生成的fpkm表达量的文件名 例如"fpkm.txt"
    :params group_file, 用户输入的样本的分组信息
    :param: venn_table_path, outfile, venn_table文件路径 当更新样本名时，输出样本的venn图分组信息表, 当样本分组数目小于6组时，默认比较样本的全部分组信息；当大于6时，默认只取前6组分组信息；默认为None
    """
    path_out = os.path.split(file_path_count)[0]
    out_path_count = os.path.join(path_out, out_count)
    out_group_file = os.path.join(path_out, "edger_group")
    if os.path.exists(file_path_count):
        pass
    else:
        raise Exception("{} 输入文件不存在".format(file_path_count))
    sample_name={}
    with open(file_path_count, "r+") as f, open(out_path_count, "w+") as ss, open(group_file, "r+") as l, open(out_group_file, "w+") as m:
        old_name =f.readline().strip().split("\t")[1:]
        #print old_name
        new_name ={}
        sample_condition=[]
        condition=[]
        for group in l:
            if group.startswith("#"):
                next
            else:
                group_line = group.strip().split("\t")
                if len(group_line)==2:
                    condition.append(group_line[1])
                    sample_condition.append(group_line[0])
        index={}
        for i in condition:
            index[i]=index.get(i,0)+1
        #print index
        for keys in range(len(index.keys())): #样本sample1-n选择
            _keys = index.keys()[keys]
            id =[]
            #print _keys
            for j in range(len(condition)): #样本同属一个类别的index
                if _keys ==condition[j]:
                    id.append(j)
            if id:
                for s in range(len(id)): #给样本设置新的名称
                    _new_name="sample"+str(keys+1)+"_"+str(s+1)
                    print _new_name
                    new_name[sample_condition[id[s]]]=_new_name
            else: 
                raise Exception("没有找到分类标签对应的样本名称 请检查")
        #print new_name
        ss.write("geneid"+"\t"+"\t".join(new_name.values())+"\n")
        m.write("#sample"+"\t"+"group"+"\n")
        for keys in new_name.keys(): # new sample name and condition
                id_type=sample_condition.index(keys)
                cond=condition[id_type]
                m.write(new_name[keys]+"\t"+cond+"\n")
        for ll in f:
            line=ll.strip().split("\t")
            ss.write("\t".join(line)+"\n")
    #if os.path.exists(file_path):
    #    os.remove(file_path)
    #    shutil.copy2(out_path, file_path)
    if os.path.exists(out_path_count):
        fpkm=change_fpkm_file(out_path_count,file_path_fpkm,out_fpkm)
    else:
        raise Exception("没有根据分类信息生成新的样本文件，请检查！")
    print "end"
    if venn_table_path:
        venn_new_path = venn_change(out_group_file, venn_table_path)
        return out_path_count, new_name, out_group_file, fpkm, venn_table_path
    else:
        return out_path_count, new_name, out_group_file, fpkm

def venn_change(out_group_file, venn_table_path):
    with open(out_group_file,"r+") as p1, open(venn_table_path, "w+") as p2:
         line1=p1.readline().strip()
         p2.write("#control\tgroup\n")
         if line1.startswith("#"):
             pass
         store_sample_info = {}  #分组对应的样本信息
         #ample_group_num=[]  #提取出所有的分组
         for line in p1:
             line2=line.strip().split("\t")
             if line2[1] not in store_sample_info.keys():
                 store_sample_info[line2[1]]={}
                 store_sample_info[line2[1]]=[line2[0]]
             else:
                 store_sample_info[line2[1]].append(line2[0])    
         if len(store_sample_info.keys())>6:  #样本分组数目大于6，则只取前6组；否则，取全部分组；
             p=0
             for keys, values in store_sample_info.items():
                 p+=1
                 if p<=6:
                    for k in values:
                        p2.write(k+"\t"+keys+"\n")
                 else:
                     break
         else:
             for keys, values in store_sample_info.items():
                  for k in values:
                       p2.write(k+"\t"+keys+"\n")
    return venn_table_path
    #return out_path_count, new_name, out_group_file, fpkm, venn_table_path

def change_fpkm_file(out_path_count,file_path_fpkm,out_fpkm):
    out_path_fpkm=os.path.join(os.path.split(file_path_fpkm)[0],out_fpkm)
    with open(out_path_count,"r+") as f,open(file_path_fpkm,"r+") as f1,open(out_path_fpkm,"w+") as f2:
        header=f.readline().strip().split("\t")
        f2.write("\t".join(header)+"\n")
        f1.readline()
        for line in f1:
            line1=line.strip().split("\t") 
            f2.write("\t".join(line1)+"\n")
        return out_path_fpkm



if __name__ == "__main__":
      a, b, c,d,e  = changesamplename(file_path_count="/mnt/ilustre/users/sanger-dev/workspace/20161130/RefrnaAssemble_zebra_single_ref_rna/Express/output/mergefeaturecounts_express/featurecounts_count.txt", file_path_fpkm="/mnt/ilustre/users/sanger-dev/workspace/20161130/RefrnaAssemble_zebra_single_ref_rna/Express/output/mergefeaturecounts_express/featurecounts_fpkm.txt", out_count="count.txt", out_fpkm="fpkm.txt", group_file="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/workflow/test_files/fastq/zebrafish/edger_group", venn_table_path = "/mnt/ilustre/users/sanger-dev/workspace/20170221/Single_htseq_sample_v2/Htseq/output/venn_table.xls")
      print c
      print d
      print "开始输出样本置换的名称信息!"
      print b
      print e
