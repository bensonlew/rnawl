# -*- coding: utf-8 -*-
# __author__ = 'zzg'

import os
import re
import shutil
from biocluster.core.exceptions import FileError
import pandas as pd
import HTMLParser
from biocluster.core.exceptions import OptionError
from biocluster.file import download,exists


def fasta_dir_sort(file_dir,sort_dir):
    #fasta_sort_dir = self.work_dir + "/fasta_dir_sort"
    if os.path.exists(sort_dir):
        shutil.rmtree(sort_dir)
    os.mkdir(sort_dir)
    if os.path.exists(file_dir + '/' + "list.txt"):
        sample_dict = {}
        with open(file_dir + '/' + "list.txt") as t:
            for i in t.readlines()[1:]:
                if i:
                    if len(i.split("\t")) < 2:
                        raise Exception("list文件应为两列数据！")
                    if i.strip().split("\t")[1].endswith(".gz"):
                        if os.path.exists(file_dir + '/' + i.strip().split("\t")[1]):
                            os.system("gzip -d {}".format(file_dir + '/' + i.strip().split("\t")[1]))
                        else:
                            raise Exception("文件{}不存在！".format(i.strip().split("\t")[1]))
                    if i.split("\t")[0] not in sample_dict.keys():
                        sample_dict[i.split("\t")[0]] = [i.strip().split("\t")[1].strip(".gz")]
                    else:
                        sample_dict[i.split("\t")[0]].append(i.split("\t")[1].strip(".gz"))
    else:
        raise Exception("list文件不存在！")
    for sample in sample_dict.keys():
        if len(sample_dict[sample]) > 1:
            cat_file = sort_dir + '/' + sample + '.fasta'
            os.system("cat {} >> {}".format(" ".join(sample_dict[sample]), cat_file))
        else:
            os.link(file_dir + '/' + sample_dict[sample][0], sort_dir + '/' + sample + '.fasta')
    return sample_dict

def fasta_dir_sort_faa(file_dir,sort_dir):
    #fasta_sort_dir = self.work_dir + "/fasta_dir_sort"
    if os.path.exists(sort_dir):
        shutil.rmtree(sort_dir)
    os.mkdir(sort_dir)
    if os.path.exists(file_dir + '/' + "list.txt"):
        sample_dict = {}
        with open(file_dir + '/' + "list.txt") as t:
            for i in t.readlines()[1:]:
                if i:
                    if len(i.split("\t")) < 2:
                        raise Exception("list文件应为两列数据！")
                    if i.strip().split("\t")[1].endswith(".gz"):
                        if os.path.exists(file_dir + '/' + i.strip().split("\t")[1]):
                            os.system("gzip -d {}".format(file_dir + '/' + i.strip().split("\t")[1]))
                        else:
                            raise Exception("文件{}不存在！".format(i.strip().split("\t")[1]))
                    if i.split("\t")[0] not in sample_dict.keys():
                        sample_dict[i.split("\t")[0]] = [i.strip().split("\t")[1].strip(".gz")]
                    else:
                        sample_dict[i.split("\t")[0]].append(i.split("\t")[1].strip(".gz"))
    else:
        raise Exception("list文件不存在！")
    for sample in sample_dict.keys():
        if len(sample_dict[sample]) > 1:
            cat_file = sort_dir + '/' + sample + '_CDS.faa'
            os.system("cat {} >> {}".format(" ".join(sample_dict[sample]), cat_file))
        else:
            os.link(file_dir + '/' + sample_dict[sample][0], sort_dir + '/' + sample + '_CDS.faa')
    return sample_dict

def gfa_dir_sort(file_dir,sort_dir):
    #fasta_sort_dir = self.work_dir + "/fasta_dir_sort"
    if os.path.exists(sort_dir):
        shutil.rmtree(sort_dir)
    os.mkdir(sort_dir)
    if os.path.exists(file_dir + '/' + "list.txt"):
        sample_dict = {}
        with open(file_dir + '/' + "list.txt") as t:
            for i in t.readlines()[1:]:
                if i:
                    if len(i.split("\t")) < 2:
                        raise Exception("list文件应为两列数据！")
                    if i.split("\t")[0] not in sample_dict.keys():
                        sample_dict[i.split("\t")[0]] = [file_dir + '/' + i.split("\t")[1].strip()]
                    else:
                        sample_dict[i.split("\t")[0]].append(file_dir + '/' + i.split("\t")[1].strip())
    else:
        raise Exception("list文件不存在！")
    for sample in sample_dict.keys():
        if len(sample_dict[sample]) > 1:
            cat_file = sort_dir + '/' + sample + '.gfa'
            os.system("cat {} >> {}".format(" ".join(sample_dict[sample]), cat_file))
        else:
            os.link(sample_dict[sample][0], sort_dir + '/' + sample + '.gfa')
    return sample_dict

# 将行名用name1代替
def rename_name(raw_file,new_file):
    name_dict = {}
    if os.path.exists(new_file):
        os.remove(new_file)
    #os.rename(file,file+"_raw_name")
    with open(raw_file) as f,open(new_file,"w") as t:
        data = f.readlines()
        t.write(data[0])
        i = 0
        for line in data[1:]:
            i+=1
            name_dict["tmpxx_name"+str(i) + "a"] = line.strip().split("\t")[0]
            t.write("tmpxx_name"+str(i) + "a" + "\t" + "\t".join(line.strip().split("\t")[1:]) + "\n")
    return name_dict

# 将行名替换回来
def rename_name_back(file,name_dict,index=0,header=True):
    os.rename(file, file + "_raw_name")
    with open(file + "_raw_name") as f,open(file,"w") as t:
        data = f.readlines()
        if header:
            t.write(data[0])
            data1 = data[1:]
        else:
            data1 = data
        for line in data1:
            if line.strip().split("\t")[index] in name_dict:
                t.write(("\t".join(line.strip().split("\t")[0:index]) + "\t" + name_dict[line.strip().split("\t")[index]] + "\t" + "\t".join(line.strip().split("\t")[index+1:])).strip() + "\n")
            elif line.strip().split("\t")[index] in name_dict.values():
                t.write(line.strip() + "\n")
            else:
                os.remove(file)
                os.rename(file + "_raw_name",file)
                raise FileError("功能{}在结果文件中未找到！".format(line.strip().split("\t")[0]))
    os.remove(file + "_raw_name")

#csa
def translate_table(self, table, out_table):
    """
    转置一个表
    :param table:输入表path
    :param out_table: 输出表path
    :return:
    """
    self.logger.info("开始对table表进行转置和修改结果表！")
    data = pd.read_table(table, sep='\t', header=0)
    name_list = data.columns
    name_list2 = ['name']
    for i in name_list[1:]:
        name_list2.append(i)
    data.columns = name_list2
    data2 = data.T
    data2 = data.astype('str').T
    data2.to_csv(out_table, sep='\t', header=0)


def meta_get_file(project_data_raw,work_dir):
    """
    meta打通小工具picrust2获取文件
    :param table:输入表path
    :param out_table: 输出表path
    :return:
    """
    analysis_type_dict = {"Enzyme": "KEGG/prediction_enzyme.xls",
                          "KO": "KEGG/prediction_KO.xls",
                          "Module": "KEGG/prediction_module.xls",
                          "Pathway_Level1": "KEGG/prediction_pathway.L1.xls",
                          "Pathway_Level2": "KEGG/prediction_pathway.L2.xls",
                          "Pathway_Level3": "KEGG/prediction_pathway.L3.xls",
                          "MetaCyc_pathway": "MetaCyc/MetaCyc_pathway_pred.xls",
                          "COG": "COG/prediction_cog.xls",
                          "Function": "COG/prediction_function.xls",
                          "funguild": "FUNGuild_guild.txt"
                          }
    input_file = work_dir + "/input_table_raw.txt"
    group_table = work_dir + "/group_table.txt"

    project_data = eval(HTMLParser.HTMLParser().unescape(project_data_raw))
    analysis_type = project_data["analysis_type"]
    del_list = []
    print(project_data["path"] + "/" + analysis_type_dict[analysis_type])
    download(project_data["path"] + "/" + analysis_type_dict[analysis_type], work_dir + "/download.txt")
    if analysis_type in ["Enzyme", "KO", "Module", "Pathway_Level2", "MetaCyc_pathway", "COG", "Function"]:
        del_list = [1]
    elif analysis_type == "Pathway_Level1":
        del_list = []
    elif analysis_type == "Pathway_Level3":
        del_list = [0, 2, 3]
    elif analysis_type == "funguild":
        pass
    else:
        raise OptionError('分析方法不正确！')

    all_sample = []
    group_data = project_data["group_dict"]
    with open(group_table, "w") as t:
        t.write("#name\tgroup\n")
        for x in group_data:
            for y in group_data[x]:
                t.write(y + "\t" + x + "\n")
                all_sample.append(y)
    if analysis_type != "funguild":
        data = pd.read_table(work_dir + "/download.txt", sep="\t")
        data.drop(data.columns[del_list], axis=1, inplace=True)
    else:
        data = pd.read_table(work_dir + "/download.txt", sep="\t")
        #for i in range(len(data.columns.values)):
        #    if data.columns.values[i] == "taxonomy":
        #        del_list = list(data.columns.values[i:])
        #data.drop(del_list, axis=1, inplace=True)
    if "abundance_top" in project_data:
        data["all_sum"] = data.apply(lambda x: x[1:].sum(), axis=1)
        data = data.sort_values(by='all_sum', ascending=False)[:int(project_data["abundance_top"])]
        data.drop("all_sum", axis=1, inplace=True)
    for i in data.columns[1:]:
        if i not in all_sample:
            data.drop(i, axis=1, inplace=True)
    data.to_csv(input_file, sep="\t", index=None)
    return input_file,group_table
