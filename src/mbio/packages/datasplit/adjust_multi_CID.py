#!/usr/bin/env python
# coding: utf-8
# 还是得对行运用函数；

import pandas as pd
import numpy as np
from collections import Counter
import os



class DealMultiCid(object):
    '''对于一个代谢物名称对应多个CID的，按照一定的规则进行过滤'''
    def __init__(self,data,cid2name,cid2map,cid2detail,out):
        self.data = data
        self.cid2name = cid2name
        self.cid2map = cid2map
        self.cid2detail = cid2detail
        self.out = out
        self.KEGG_list = pd.read_csv(self.cid2name,sep="\t",names=["CID","NAME"], header=None)
        self.KEGG_list.loc[:,"CID"] = self.KEGG_list.loc[:,"CID"].str.replace("cpd:","")
        self.long_KEGG = self.KEGG_list.drop('NAME', axis=1).join(self.KEGG_list.loc[:,'NAME'].str.split('; ', expand=True).stack().reset_index(level=1, drop=True).rename('NAME'))
        self.L_STRUCT_list = list(set(self.long_KEGG[self.long_KEGG.loc[:,"NAME"].str.contains('^L-')]["CID"].tolist()))

        self.KEGG_MASS = pd.read_csv(self.cid2detail,sep="\t",usecols=["ENTRY","EXACT_MASS"])
        self.KEGG_MASS.columns = ["CID","EXACT_MASS"]
        self.CPD_PATH = pd.read_csv(self.cid2map,sep="\t",header=None,names=["CID","PATHS"])

        self.CPD_PATH.loc[:,'CID'] = self.CPD_PATH.loc[:,'CID'].str.replace("cpd:",'')
        self.CPD_2_PATH = self.CPD_PATH.groupby(by ="CID").agg(';'.join)
        self.CPD_2_PATH = self.CPD_2_PATH.reset_index()

    def deal_each_meta(self):
        data_df = pd.read_csv(self.data,sep="\t",quoting=0, quotechar='"', encoding='utf-8')
        mt_data = data_df[data_df.loc[:,"KEGG Compound ID"].str.contains(";")]
        other_data = data_df[~data_df.loc[:,"KEGG Compound ID"].str.contains(";")]

        extend_data = mt_data.drop('KEGG Compound ID', axis=1).join(mt_data.loc[:,'KEGG Compound ID'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('KEGG Compound ID'))
        extend_data_mass = pd.merge(left=extend_data,right=self.KEGG_MASS,how="left",left_on="KEGG Compound ID",right_on="CID",sort=False).drop(["CID"],axis=1)
        extend_data_mass_path = pd.merge(left=extend_data_mass,right=self.CPD_2_PATH,how="left",left_on="KEGG Compound ID",right_on="CID",sort=False).drop(["CID"],axis=1)
        extend_data_mass_path_names = pd.merge(left=extend_data_mass_path,right=self.KEGG_list,how="left",left_on="KEGG Compound ID",right_on="CID",sort=False).drop(["CID"],axis=1)
        extend_data_mass_path_names.loc[:,"configuration"] = np.where(extend_data_mass_path_names["KEGG Compound ID"].isin(self.L_STRUCT_list), 'L_STRUCT', 'Non_L')

        def judge_mass_cal(row):
            row_mass = row[["m/z","EXACT_MASS"]]
            if row_mass.shape[0] <2:
                row_mass = row[['mass',"EXACT_MASS"]]
                row_mass = row_mass.astype(float)
                delta = abs(row_mass["mass"] - row_mass["EXACT_MASS"])
            else:
                row_mass = row_mass.astype(float)
                delta = abs(row_mass["m/z"] - row_mass["EXACT_MASS"])
            return delta

        def judge_group(group):
            group.loc[:,"delta"] = group.apply(judge_mass_cal,axis=1)
            group_sort = group.sort_values(by='delta',ascending=True)
            group_delta = group_sort["delta"].tolist()
            min_num = min(group_delta)
            group_times = Counter(group_delta)
            other_times = len(group_delta) - group_times[min_num]
            group_sort.loc[:,"Condition"] = group_times[min_num]*["KEEP"] + other_times * ["DELETE"]
            group_end = group_sort.loc[group_sort.loc[:,"Condition"].isin(["KEEP"]),:]
            return group_end

        def judge_name(row):
            result = "DELETE"
            mj_name = row["Metabolite"]
            kegg_names = row["NAME"].split('; ')
            lower_names = [x.lower() for x in kegg_names]
            if mj_name.lower() in lower_names:
                result = 'KEEP'
            return result

        def judge_path(row):
            result ="DELETE"
            tmp_df =row["PATHS"]
            if tmp_df is np.nan:
                result ="DELETE"
            else:
                result="KEEP"
            return result

        def judge_struct(row):
            result = "DELETE"
            if row["configuration"] == "L_STRUCT":
                result = "KEEP"
            return result

        def judge_other3(group):
            group.loc[:,"Condition_PA"] = group.apply(judge_path,axis=1)
            group.loc[:,"Condition_ST"] = group.apply(judge_struct,axis=1)
            group.loc[:,"Condition_NM"] = group.apply(judge_name,axis=1)
            pa_group = group.loc[group.loc[:,"Condition_PA"].isin(["KEEP"]),:]
            if pa_group.shape[0] ==1:
                pa_group.loc[:,"Condition"] = "KEEP"
                end_group= pa_group.drop(["Condition_ST","Condition_PA","Condition_NM"],axis=1)
                return end_group
            else:
                st_group = group.loc[group.loc[:,"Condition_ST"].isin(["KEEP"]),:]
                if st_group.shape[0] ==1:
                    st_group.loc[:,"Condition"] = "KEEP"
                    end_group= st_group.drop(["Condition_ST","Condition_PA","Condition_NM"],axis=1)
                    return end_group
                else:
                    nm_group = group.loc[group.loc[:,"Condition_NM"].isin(["KEEP"]),:]
                    if nm_group.shape[0] ==1:
                        nm_group.loc[:,"Condition"] = "KEEP"
                        end_group= nm_group.drop(["Condition_ST","Condition_PA","Condition_NM"],axis=1)
                        return end_group
                    else:
                        tmp_group = group.drop(["Condition_ST","Condition_PA","Condition_NM"],axis=1)
                        end_group = tmp_group.iloc[0:1,:]
                        end_group.loc[:,"KEGG Compound ID"] = ';'.join(tmp_group.loc[:,"KEGG Compound ID"].tolist())
                        end_group.loc[:,"Condition"] = "KEEP"
                        return end_group

        s1_data = extend_data_mass_path_names.groupby(by=["Metabolite"])
        end_results = []
        for name,group in s1_data:
            group.loc[:,"EXACT_MASS"]= group.loc[:,"EXACT_MASS"].replace('None',0)
            st_1_group = judge_group(group)
            if st_1_group.shape[0] >1:
                st_end_group = judge_other3(st_1_group)
                if st_end_group.shape[0] == 0:
                    print(group)
                    break
                end_results.append(st_end_group)
            else:
                end_results.append(st_1_group)
        # print(end_results)
        end_df = pd.concat(end_results, axis=0, sort=False)
        end_mt = end_df.drop(["Condition","configuration","PATHS","NAME","EXACT_MASS","delta"], axis =1)
        # print("----------")
        # print(other_data)
        # print(end_mt)
        end_data = pd.concat([other_data,end_mt],axis =0,sort=False)
        end_data.to_csv(self.out, sep="\t", quoting=0, quotechar='"', index=False, header=True, encoding='utf-8')
        end_data.to_excel(self.out + 'x', index=False, header=True,encoding='utf-8')


def adjust_mCID(file, out):
    # cid2name = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\CID2NAME.20210409.v98.0.txt'
    # cid2map = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\Compound_Pathway.list'
    # cid2detail = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\KEGG_compound.20210526_addHMDB.xls'
    # print(F'开始处理{file}里的多个cid的情况')
    print('开始处理%s里的多个cid的情况' % file)
    need_deal = DealMultiCid(file, cid2name, cid2map, cid2detail, out)
    try:
        need_deal.deal_each_meta()
    except Exception:
        # print(F'{file}不需要处理')
        print('%s不需要处理' % file)
        os.rename(file, out)


if __name__ == '__main__':
    # cid2name = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\CID2NAME.20210409.v98.0.txt'
    # cid2map = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\Compound_Pathway.list'
    # cid2detail = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\KEGG_compound.20210526_addHMDB.xls'
    need_deal = DealMultiCid(args.data, args.cid2name, args.cid2map, args.cid2detail,args.out)
    need_deal.deal_each_meta()
