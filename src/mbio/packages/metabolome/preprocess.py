# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import numpy as np
import math
import pandas as pd
from rfp import missrf
import copy

class Preprocess(object):

    def __init__(self, file, index_col, control_sample=None,inner_file=None,group_file=None):
        self.data_in_treat = pd.read_table(file, index_col=index_col)
        self.control_sample = control_sample
        #self.stat_data = ""  # metab_id rsd  sum
        self.rsd = None
        if group_file:
            self.group_map = self.get_group_map(group_file)
        else:
            self.group_map = None
        if inner_file:
            self.inner_data = pd.read_table(inner_file, index_col=index_col)
            self.has_inner_data = True
        else:
            self.has_inner_data = False


    def get_group_map(self,group_file):
        tmp_dic = {}
        with open(group_file) as fr:
            fr.readline()
            for line in fr:
                line = line.strip()
                if line =='':
                    continue
                spline = line.split('\t')
                if spline[1] not in tmp_dic:
                    tmp_dic[spline[1]] = []
                tmp_dic[spline[1]].append(spline[0])
        return tmp_dic




    def fillna(self, method, na=0):
        if method == "rf":
            self.data_in_treat = self.data_in_treat.apply(missrf, axis=1, args=(self.data_in_treat,))
            return
        # 随机森林之外的处理方法
        for index in self.data_in_treat.index:
            if self.data_in_treat.loc[index].max() == 0 and na == 0:
                raise Exception("%s行值均为0" % index)
            value = self.get_value(self.data_in_treat.loc[index], method)
            self.data_in_treat.loc[index] = self.data_in_treat.loc[index].replace(na, value)

    def get_value(self, ser, method, na=0):
        value = na
        ser = ser.replace(value, np.nan)
        if method == 'min':
            value = ser.min()
        elif method == 'max':
            value = ser.max()
        elif method == 'mean':
            value = ser.mean()
        elif method == 'median':
            value = ser.median()
        elif method == 'rf':
            value = ser.sample().values[0]
            while value == na or np.isnan(value):
                value = ser.sample().values[0]
        elif method == 'none':
            pass
        return value

    def control(self, method):
        if len(self.control_sample) < 2:
            raise Exception("质控样本不能少于2")
        tmp_data = self.data_in_treat[self.control_sample].astype(float)
        self.data_in_treat['rsd'] = tmp_data.apply(lambda x: x.std()/x.mean(), axis=1)
        self.rsd = self.data_in_treat['rsd']
        self.has_rsd = True
        self.data_in_treat['sum'] = tmp_data.apply(lambda x: x.sum(), axis=1)
        self.data_in_treat = self.data_in_treat[self.data_in_treat['rsd'] <= method]
        if len(self.data_in_treat) == 0:
            raise Exception("质控后代谢物数量不足")
        remain_index = self.data_in_treat.index
        return remain_index

    def norm(self, method):
        def norm_value(df, scale):
            tmp_df = df
            if "rsd" in df.index:
                tmp_df = df[:-2]
            for index, i in enumerate(tmp_df):
                i = i/scale[index]
                df[index] = i
            return  df
        scale = 0
        if method in ["sum", "median", "mean"]:
            scale = self.norm_stat(method)
        elif "sample" in method:
            method = method.lstrip("sample:")
            scale = self.norm_sample(method)
        elif "inner" in method:
            method = method.lstrip("inner:")
            scale = self.norm_inner(method)
        elif method == "none":
            return
        self.data_in_treat = self.data_in_treat.apply(norm_value, axis=1, args=(scale,))

    def norm_inner(self, inner):
        if self.has_inner_data:
            inner_value = self.inner_data.loc[inner]
        else:
            inner_value = self.data_in_treat.loc[inner]
        if 0 in inner_value.tolist():
            raise Exception("内参代谢物%s含有0" % inner)
        for header in ["rsd", "sum"]:
            if header in inner_value:
                inner_value = inner_value.drop([header])
        ref_value = inner_value.max()
        scale = inner_value / ref_value
        return scale

    def norm_sample(self, sample):
        sum = self.data_in_treat.apply(lambda x: x.sum())
        ref_value = sum[sample]
        scale = sum / ref_value
        return scale

    def norm_stat(self, method):
        static = ""
        if method == "mean":
            static = self.data_in_treat.apply(lambda x: x.mean())
        elif method == "median":
            static = self.data_in_treat.apply(lambda x: x.median())
        elif method == "sum":
            static = self.data_in_treat.apply(lambda x: x.sum())
        for header in ["rsd", "sum"]:
            if header in static:
                static = static.drop([header])
        ref_value = static.max()
        scale = static / ref_value
        return scale

    def log(self, method):
        # 判断表中是否存在rsd和sum两列，如果有去除后再进行log处理
        self.column_list = self.data_in_treat.columns.tolist()
        if "rsd" in self.column_list and "sum" in self.column_list:
            self.dataframe = self.data_in_treat.drop(["rsd","sum"], axis=1)
        else:
            self.dataframe = self.data_in_treat
        if method == 2:
            #self.data_in_treat = np.log2(self.data_in_treat + 1)
            data_in_log = np.log2(self.dataframe + 1)
            if "rsd" in self.column_list and "sum" in self.column_list:
                self.data_in_treat = pd.concat([data_in_log, self.data_in_treat["rsd"], self.data_in_treat["sum"]], axis=1)
            else:
                self.data_in_treat = data_in_log
        elif method == 10:
            #self.data_in_treat = np.log10(self.data_in_treat + 1)
            data_in_log = np.log10(self.dataframe + 1)
            if "rsd" in self.column_list and "sum" in self.column_list:
                self.data_in_treat = pd.concat([data_in_log, self.data_in_treat["rsd"], self.data_in_treat["sum"]], axis=1)
            else:
                self.data_in_treat = data_in_log
        elif method == "none":
            return
        else:
            self.log_num(method)

    def log_num(self, num):
        def get_log(df, number):
            for index,i in enumerate(df):
                value = math.log(i, number)
                df[index] = value
            return df
        new_data = (self.dataframe + 1).apply(get_log, axis=1, args=(num,))
        if "rsd" in self.column_list and "sum" in self.column_list:
            self.data_in_treat = pd.concat([new_data, self.data_in_treat["rsd"], self.data_in_treat["sum"]], axis=1)
        else:
            self.data_in_treat = new_data
        #new_data = (self.data_in_treat + 1).apply(get_log, axis=1, args=(num,))
        #self.data_in_treat = new_data
        # self.data_in_treat = self.data_in_treat.apply(lambda x: math.log(x, num))

    def to_table(self, outfile):
        self.data_in_treat.to_csv(outfile, sep="\t", quoting=3)

    def pipline(self, piplist, pipdic):
        """
        整体的处理流程
        :param piplist: ["fillna", "control", "norm", "log", "to_table"]
        :param pipdic: {"fillna": "mean", "control": "0.3", "rm_nan":0.5,"rm_type":"all"}
        :return:
        """
        remain_index = self.data_in_treat.index
        for ele in piplist:
            if ele == "fillna":
                self.fillna(pipdic[ele])
            elif ele == "rm_nan":  ##zouguanqing 20190604
                #remain_index = self.rm_nan(pipdic[ele])
                if 'rm_type' not in pipdic:
                    remain_index = self.rm_nan(pipdic[ele]) #v3
                else:
                    if pipdic["rm_type"] == 'each':
                        remain_index = self.rm_nan_by_group(pipdic[ele]) #v3
                    else:
                        remain_index = self.rm_nan(pipdic[ele]) #v3

            elif ele == "norm":
                self.norm(pipdic[ele])
            elif ele == "control" and pipdic['control'] != 'none':
                remain_index = self.control(pipdic[ele])
            elif ele == "log":
                if 'control' not in piplist:
                    self.get_rsd()
                self.log(pipdic[ele])
            elif ele == "to_table":
                if 'log' not in piplist:
                    self.get_rsd()
                self.to_table(pipdic[ele])
                if self.has_rsd:
                    self.rsd.to_csv(pipdic[ele]+'_rsd.txt',sep='\t')
        return remain_index

    def get_rsd(self):
        if len(self.control_sample) < 2:
            self.rsd =None
            self.has_rsd = False
            return False
            #raise Exception("质控样本不能少于2")
        tmp_data = self.data_in_treat[self.control_sample].astype(float)
        tmp_data['rsd'] = tmp_data.apply(lambda x: x.std()/x.mean(), axis=1)
        self.rsd = tmp_data['rsd']
        self.has_rsd = True

    def rm_nan_by_group(self, threshold):
        if not self.group_map:
            raise Exception('去除缺失值当选择 每一组时，要提供分组文件')

        rm_qc_group = copy.deepcopy(self.group_map)
        if 'QC' in rm_qc_group:
            rm_qc_group.pop('QC')
        def cent_min_group_fun(one):
            min_loss_per = 1
            for g in rm_qc_group:
                if g == 'QC':
                    continue
                loss = 0
                all = 0
                for s in rm_qc_group[g]:
                    all +=1
                    if one[s]==0 or one[s]=='-':
                        loss +=1
                loss_per = float(loss)/all
                if loss_per < min_loss_per:
                    min_loss_per = loss_per

            return min_loss_per

        self.data_in_treat = self.data_in_treat[self.data_in_treat.apply(lambda x: cent_min_group_fun(x) <= threshold,axis=1)]
        return self.data_in_treat.index



    def rm_nan(self,threshold):  ## 20190604
        def cent_fun(values):
            n = 0
            all = 0
            for v in values:
                if v==0 or v=='-':
                    n+=1
                all +=1
            return float(n)/all

        if self.control_sample:  ## 去掉缺失值不包括qc样本
            all_columns = self.data_in_treat.columns.tolist()
            rm_control = []
            for i in all_columns:
                if i not in self.control_sample:
                    rm_control.append(i)
            self.rm_control_data = self.data_in_treat[rm_control]
            self.data_in_treat = self.data_in_treat[self.rm_control_data.apply(lambda x: cent_fun(x) <= threshold,axis=1)]

        else:
            self.data_in_treat = self.data_in_treat[self.data_in_treat.apply(lambda x: cent_fun(x) <= threshold,axis=1)]
        return self.data_in_treat.index
