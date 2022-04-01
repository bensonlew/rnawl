# -*- coding: utf-8 -*-
# __author__ = "shaohua.yuan"
# last_modifiy = modified 2018.06.13


import os, re, sys, argparse, glob
import pandas as pd
import re
import copy

def merge_table(table1, table2, merge_columns, outfile=None, how="inner", out_columns=[]):
    """
    根据table中的某一列合并两个table文件
    :param table1: 输入的table文件
    :param table2: 输入的table文件
    :param merge_columns: 合并的columns
    :param how: inner 交集，outer 并集，left 取左侧的行，right 取右侧的行
    :param out_columns：
    :return: 返回生成的合并文件outfile
    """
    if isinstance(table1, str):
        t1 = pd.read_table(table1, sep="\t", header=0)
    else:
        t1 = table1
    if isinstance(table2, str):
        t2 = pd.read_table(table2, sep="\t", header=0)
    else:
        t2 = table2
    merge_table = pd.merge(t1, t2, how=how, on=merge_columns)
    if out_columns:
        merge_tale = merge_table[out_columns]
    merge_table = merge_table.fillna("-")
    if outfile:
        merge_table.to_csv(outfile, sep="\t", index=False, quoting=3)
    return merge_table

def get_sub_table(table, sub_columns, outfile=None,has_name_set=None):
    #df = pd.read_table(table, sep="\t", header=0)
    df = pd.DataFrame(table)
    if has_name_set:
        df = df[df[sub_columns].apply(lambda x: x in has_name_set)]  # 筛选有代谢物名称的metab id
    df = pd.DataFrame(df[sub_columns])

    if outfile:
        df.to_csv(outfile, sep="\t", index=False)
    return df




def  _filter_talbe(table,name,type,value):
    if type == 'gt':
        filter = table[table[name] > value]
    elif type == 'lt':
        filter = table[table[name] < value]
    elif type == 'gte':
        filter = table[table[name] >= value]
    elif type == 'lte':
        filter = table[table[name] <= value]
    elif type == 'eq':
        filter = table[table[name] == value]
    else:
        return table
    return filter
    #filter_conditons = [
    # {'P_value',{'type':'>','value':0.05},},  #key : P_value, Vip_oplsda, Vip_plsda, fdr FC

def filter_table(table, outfile=None,filter_conditions=None):  #zouguanqing
    if filter_conditions:
        filter1 = ''
        filter2 = ''
        if 'FC_up' in filter_conditions.keys():
            condition = filter_conditions['FC_up']
            type = condition['type']
            value = condition['value']
            filter1 = _filter_talbe(table,'FC',type,value)
        if 'FC_down' in filter_conditions.keys():
            condition = filter_conditions['FC_down']
            type = condition['type']
            #value_new = 1/float(condition['value'])
            type_map = {'gt':'lt','lt':'gt','gte':'lte','lte':'gte'}
            #type_new = type_map[type]
            value_new = condition['value']
            type_new = type
            filter2 = _filter_talbe(table,'FC',type_new,value_new)

        if not isinstance(filter1,str)  and not isinstance(filter2,str):
            filter = pd.concat([filter1,filter2],axis=0)
        elif not isinstance(filter2,str):
            filter = filter2
        elif not isinstance(filter1,str):
            filter = filter1
        else:
            filter = table

        for name in filter_conditions.keys():
            if name in ['FC_up','FC_down']:
                continue
            condition = filter_conditions[name]
            type = condition['type']
            value = condition['value']
            filter = _filter_talbe(filter,name,type,value)

    else:
        filter = table[table["P_value"] < 0.05]
        filter = filter[filter["Vip_oplsda"] > 1]
    return filter

def run_merge(test_dir, pls_dir, name, outDir):
    all_test_file = os.listdir(test_dir)
    all_file_dir = os.listdir(pls_dir)
    groups = name.split(";")
    group_merge_file = {}
    for each in groups:
        group = each
        test_file = glob.glob(test_dir + "/" + group + "_result.xls")
        test_table = get_table(test_file[0])
        plsda_table = get_table(pls_dir + "/" + group + "/PLS-DA.vip.xls")
        oplsda_table = get_table(pls_dir + "/" + group + "/OPLS-DA.vip.xls")
        test_table.rename(columns={test_table.columns[0]: "Metab"}, inplace=True)
        plsda_table.columns = ["Metab", "Vip_plsda"]
        oplsda_table.columns = ["Metab", "Vip_oplsda"]
        vip = merge_table(plsda_table, oplsda_table, "Metab")
        print vip.head()
        final_merge = merge_table(vip, test_table, "Metab", how="inner")
        print final_merge.head()
        if not os.path.exists(outDir + "/tmp_DiffStat"):
            os.mkdir(outDir + "/tmp_DiffStat")
        outfile = outDir + "/tmp_DiffStat/" + group + ".diff.exp.xls"
        final_merge.to_csv(outfile, sep="\t", index=False)
        group_merge_file[group] = outfile
    return group_merge_file

def creat_sub_metab(merge_file_dict, outDir,desc=None,filter_conditions=None): #desc 文件，用来提取用代谢物名称的metab id
    if desc:
        desc_data = pd.read_table(desc,sep='\t',header=0)
        has_name_metabs = desc_data[desc_data['Metabolite'].apply(lambda x: False if re.match('^metab_\d*',x) or re.match('^pos_\d*$',x) or re.match('^neg_\d*',x) else True)]['metab_id']
        has_name_set = set(has_name_metabs.tolist())

    if not os.path.exists(outDir + "/Metabset"):
        os.mkdir(outDir + "/Metabset")
    mul_metabsetlist = outDir + "/Metabset/mul.metabset.list.xls"
    with open(mul_metabsetlist, "w") as f:
        for each in merge_file_dict.keys():
            eachfile = merge_file_dict[each]
            eachtable = get_table(eachfile)
            filter = filter_table(eachtable,filter_conditions=filter_conditions)
            if not os.path.exists(outDir + "/Metabset"):
                os.mkdir(outDir + "/Metabset")
            outfile = outDir + "/Metabset/" + each + ".metabset.xls"
            if len(filter) >0:
                if desc:
                    submetab = get_sub_table(filter, "Metab", outfile,has_name_set=has_name_set)
                else:
                    submetab = get_sub_table(filter, "Metab", outfile)
                final_metabs = submetab["Metab"].tolist()
                if len(final_metabs) > 0:
                    each_list = ",".join(final_metabs)
                    f.write(each + "\t" + each_list + "\n")
    return mul_metabsetlist

def creat_common_diff_metab(mul_metabsetlist, outDir):
    common_list = []
    with open(mul_metabsetlist, "r") as infile:
        for line in infile:
            line = line.strip().split("\t")
            #metab_name =  line[0]
            metab_list = line[1].split(",")
            common_list = common_list + metab_list
    common_list = list(set(common_list))
    common_table = pd.DataFrame(common_list, columns=["metab_id"])
    if not os.path.exists(outDir + "/Metabset"):
        os.mkdir(outDir + "/Metabset")
    outfile = outDir + "/Metabset/Diff_intersection.metabset.xls"
    common_table.to_csv(outfile, sep="\t", index=False)

def get_table(file):
    table = pd.read_table(file, sep="\t", header=0)
    return table

def merge_mul_metabset(mul_file1, outfile, mul_file2=None, exp_des_file=None):
    mul_dict = {}
    if exp_des_file:
        map_dict = {}
        with open(exp_des_file, "r") as f:
            f.next()
            for line in f:
                line = line.strip().split("\t")
                metab_id = line[0]
                metab_des = line[1]
                if not map_dict.has_key(metab_id):
                    map_dict[metab_id] = metab_des
    with open(mul_file1, "r") as f1:
        for line in f1:
            line = line.strip().split("\t")
            group = line[0]
            metabs = line[1]
            if not mul_dict.has_key(group):
                metab_list = metabs.split(",")
                mul_dict[group] = metab_list
    if mul_file2:
        with open(mul_file2, "r") as f2:
            for line in f2:
                line = line.strip().split("\t")
                group = line[0]
                metabs = line[1]
                new_list = metabs.split(",")
                if mul_dict.has_key(group):
                    mul_dict[group] = mul_dict[group] + new_list
                else:
                    mul_dict[group] = new_list
    with open(outfile, "w") as f3:
        for eachgroup in mul_dict.keys():
            metab_lists = list(set(mul_dict[eachgroup]))
            if exp_des_file:
                metab_new_list = copy.copy(metab_lists)
                for each in metab_lists:
                    name = map_dict[each]
                    lower_name = name.lower()
                    if re.match('^pos_\d+$',lower_name) or re.match('^neg_\d+$',lower_name):
                    #if lower_name.startswith('pos') or lower_name.startswith('neg'):
                        metab_new_list.remove(each)
                metabs_str = ",".join(metab_new_list)
            else:
                metabs_str = ",".join(metab_lists)
            f3.write(eachgroup + "\t" + metabs_str + "\n")

def create_filter_conditions(filter_k,filter_t,filter_v):
    filter_ks = filter_k.split(',')
    filter_ts = filter_t.split(',')
    filter_vs = filter_v.split(',')
    conditions = {}
    for i in range(len(filter_ks)):
        name = filter_ks[i]
        t = filter_ts[i]
        v = filter_vs[i]
        tmp = {'name':name, 'type':t,'value':float(v)}
        conditions[name] = tmp
    return conditions



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d1', type=str, metavar="table1", required=True, help="merge table1")
    parser.add_argument('-d2', type=str, metavar="table2", required=True, help="merge table2")
    parser.add_argument('-o', type=str, metavar="outDir", required=True, help="output dir")
    parser.add_argument('--ms', action='store_true', help="creat metabset,True or False")
    parser.add_argument('-name', type=str, metavar="merge group name", required=True, help="merge file name")
    parser.add_argument('-desc', type=str, metavar="table3", required=False, help="metab desc")
    parser.add_argument('-filter_k',type=str,required=False,help='P_value, Vip_oplsda, Vip_plsda, fdr FC')
    parser.add_argument('-filter_t',type=str,required=False,help='gt, gte, lte,lt,eq')
    parser.add_argument('-filter_v',type=str,required=False,help='filter value')


    args = parser.parse_args()
    test_dir = args.d1
    pls_dir = args.d2
    output = args.o
    name = args.name
    if args.filter_k:
        conditions = create_filter_conditions(args.filter_k,args.filter_t,args.filter_v)
    else:
        conditions = None
    group_merge_dict = run_merge(test_dir, pls_dir, name, output)
    if args.ms:
        if args.desc:
            mul_metabsetlist = creat_sub_metab(group_merge_dict, output,desc=args.desc,filter_conditions=conditions)
        else:
            mul_metabsetlist = creat_sub_metab(group_merge_dict, output,filter_conditions=conditions)
        creat_common_diff_metab(mul_metabsetlist, output)  # -*- coding: utf-8 -*-
