#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 2021/5/10 14:55
# @Author  : U make me wanna surrender my soul
import matplotlib
matplotlib.use('Agg')
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", help="输入原始数据", type=str, required=True)
parser.add_argument("-g", help="输入分组文件", type=str, required=True)
parser.add_argument("-o", help="输入图片名称", type=str, required=True)
args = parser.parse_args()


# 加载数据集 数据中不能有字符串
inputfile = args.i
data = pd.read_csv(inputfile, index_col=0, sep='\t')
group_file = args.g
group_data = pd.read_csv(group_file, sep='\t')
group_target = np.array(group_data['group'].tolist())

# 使用TSNE进行降维处理。从4维降至2维。
tsne = TSNE(n_components=2, learning_rate=100).fit_transform(data)

# 设置画布的大小
plt.figure(figsize=(12, 6))
plt.scatter(tsne[:, 0], tsne[:, 1], c=group_target)
plt.colorbar()
# plt.show()
plt.savefig(args.o)
