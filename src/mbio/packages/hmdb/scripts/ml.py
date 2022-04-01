# -*- coding: utf-8 -*-
# __author__ = 'ghd'
"""
def iris_type(s):
    it = {'Iris-setosa': 0, 'Iris-versicolor': 1, 'Iris-virginica': 2}
    return it[s]


def desease_type(s):
    it = {'health': 0, 'disease': 1, 'u_c': 2}
    if s in it.keys():
        return it[s]
    else:
        return s
"""

def get_pk():
    from sklearn.externals.joblib import load
    with open("y_sav", "rb") as f1:
        y = load(f1)
    with open("prob_sav", "rb") as f2:
        prob = load(f2)
    return y,prob

def get_prop(prob,y):
    import numpy as np
    health_list = np.array([])
    desease_list = np.array([])
    for index,i in enumerate(y):
        if i == 0:
            if len(health_list) == 0:
                health_list = np.array([prob[index]])
            else:
                health_list = np.append(health_list, np.array([prob[index]]),axis=0)
        else:
            if len(desease_list) == 0:
                desease_list = np.array([prob[index]])
            else:
                desease_list = np.append(desease_list, np.array([prob[index]]),axis=0)
    return health_list,desease_list

def get_densi(data):
    from scipy.stats import gaussian_kde
    import numpy as np
    density = gaussian_kde(data)
    xs = np.linspace(0, 1, 200)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    return xs,density(xs)

def combine(prob, y, file):
    import pandas as pd
    hl, dl = get_prop(prob, y)
    hx,hd = get_densi(hl[:,1])
    dx,dd = get_densi(dl[:,1])
    df = pd.DataFrame({"x":hx, "y1": hd, "y2": dd})
    df.to_csv(file, sep="\t", index=False)


class MlGroup(object):
    def __init__(self, group_table):
        self.map_dic = {}
        self.desease_array = []
        self.desease_signal = {}
        self.group_table = group_table
        self.read_group()
        self.set_desease_signal()

    def read_group(self):
        desease_set = set()
        tmp_list = []
        with open(self.group_table, 'r') as f:
            for line in f.readlines():
                if "#" in line:
                    continue
                sample, group = line.strip().split('\t')
                self.map_dic[sample] = group
                # desease_set.add(group)
                #if group in ['Healthy', 'Cancer']:
                #    tmp_list = ['Healthy', 'Cancer']
                #else:
                #    desease_set.add(group)
                if group in ['Healthy']:
                    tmp_list = ['Healthy']
                else:
                    desease_set.add(group)
        self.desease_array = tmp_list + list(desease_set)

    def set_desease_signal(self):
        """
        desease_arrya:数组
        """
        for index, desease in enumerate(self.desease_array):
            self.desease_signal[desease] = index

    def desease_type(self, s):
        # self.set_desease_signal()
        if s in self.map_dic.keys():
            # map_s = self.map_dic[s]
            return self.desease_signal[self.map_dic[s]]
        else:
            return s

    @property
    def desease_hash(self):
        sample2id = {}
        for s in self.map_dic.keys():
            sample2id[s] = self.desease_signal[self.map_dic[s]]
        return sample2id

    def get_sample_list(self, label_list):
        sample_list = []
        for s in self.map_dic.keys():
            if self.map_dic[s] in label_list:
            # if self.map_dic[s] in ['Healthy', 'Cancer']:
                sample_list.append(s)
        return sample_list
