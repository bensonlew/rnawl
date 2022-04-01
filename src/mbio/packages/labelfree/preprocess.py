# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import numpy as np
import math
import pandas as pd
from collections import defaultdict
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score

class Preprocess(object):

    def __init__(self, exp_file, group_file, cutoffs='', method='rf', out=None):
        self.data_in_treat = pd.read_table(exp_file, index_col=0, sep='\t')
        df_group = pd.read_table(group_file, header=0, sep="\t")
        sample2group = zip(df_group['#sample'], df_group['group'])
        self.group2sample = defaultdict(list)
        for i in sample2group:
            self.group2sample[i[1]].append(i[0])
        self.g2cutoff = dict()
        cutoffs = cutoffs.split(';')
        groups = df_group['group']
        try:
            self.g2cutoff = {x[1]: float(x[2]) for x in zip(groups, cutoffs)}
        except:
            self.g2cutoff = {g: float(len(s)) / 2 for g, s in self.group2sample.items()}
        for g, s in self.group2sample.items():
            if self.g2cutoff[g] > len(s):
                self.g2cutoff[g] = float(len(s)) / 2
        self.method = method
        if out:
            self.out = out
        else:
            self.out = exp_file

    def fillzero(self):
        method = self.method
        for g, ss in self.group2sample.items():
            for ind in self.data_in_treat.index:
                ser = self.data_in_treat.loc[ind, ss]
                cutoff = self.g2cutoff[g]
                num_0 = ser.values.tolist().count(0)
                if not num_0 or num_0 > cutoff:
                    continue
                if method == 'rf':
                    pred = self.get_rf_value(ser, self.data_in_treat[ss])
                    for n, s in enumerate(ss):
                        if self.data_in_treat.loc[ind, s] == 0:
                            self.data_in_treat.loc[ind, s] = pred[n]
                else:
                    value = self.get_value(ser, method)
                    for s in ss:
                        if self.data_in_treat.loc[ind, s] == 0:
                            self.data_in_treat.loc[ind, s] = value

    def to_out(self):
        self.data_in_treat.to_csv(self.out, sep='\t', index=True, header=True)

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
        elif method == 'none':
            pass
        return value

    def get_rf_value(self, ser, df):
        def rf_pred(df, df_train):
            clf = RandomForestRegressor()
            fit = clf.fit(df_train.values.T, df.values)
            pred = clf.predict(df_train.values.T)
            prob = r2_score(df.values, pred)
            return pred, prob

        def rf_permu(df, df_train, count_size=1):
            perd, prob = rf_pred(df, df_train)
            if prob > 0.9 or count_size >= 10:
                return perd
            else:
                df = rf_permu(df, df_train, count_size + 1)
            return perd
        value = self.get_value(ser, "median")
        ser = ser.replace(0, value)
        index_list = df.index.tolist()
        index_list.remove(ser.name)
        df_train = df[df.index.isin(index_list)]
        perd = rf_permu(ser, df_train)
        return perd


if __name__ == '__main__':
    expfile = '/mnt/ilustre/users/sanger-dev/workspace/20190723/Labelfree_tsg_34919/treat_ref'
    groupfile = '/mnt/ilustre/users/sanger-dev/workspace/20190723/Labelfree_tsg_34919/remote_input/protein_group/group.txt'
    # method = 'rf'
    method = 'min'
    out = 'treat_min.txt'
    process = Preprocess(expfile, groupfile, method=method, out=out)
    process.fillzero()
    process.to_out()