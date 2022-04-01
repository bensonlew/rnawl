#!/usr/bin/env python
# coding: utf8

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import argparse
import os
import math

def Normaldistest(input_data_file,out_file):
	df_txt = pd.read_csv(input_data_file, sep = '\t')
	df_txt.columns = ["title","value"]
	#df = pd.read_table(input_data_file, sep =  '\n',names = ['value'])
	mean = df_txt['value'].mean()
	std = df_txt['value'].std()
	df_txt.sort_values(by = 'value', inplace = True)  # 排序
	df_txt_r = df_txt.reset_index(drop = False)  # 重新排序后，更新index
	n = len(df_txt_r)
	if 3 <= n <= 50:
		bin_num = n //5 + 2
	elif n <= 1000 and n > 50:
		bin_num = 30
	elif n > 1000:
		bin_num = 50
	else:
		raise OptionError('至少输入3个样本！')
	#ks检验
	D1,P1 = stats.kstest(df_txt_r['value'], 'norm', (mean, std)) 
	# Shapiro-Wilk test
	D2,P2 = stats.shapiro(df_txt_r['value'])
	skew = float(df_txt_r['value'].skew())
	kurt = float(df_txt_r['value'].kurt())
	data_output = open(out_file,'w')

	# 输出直方图数据
	histogram_data = np.histogram(df_txt_r['value'],bins = bin_num)
	Y1 = histogram_data[0]
	X1 = histogram_data[1]
	for i in range(len(X1)-1):
		X1[i] = (X1[i]+X1[i+1])/2
	X1 = X1[:(len(X1)-1)]
	a = 0
	b = 0
	for i in range(len(Y1)):
		a += Y1[i]
	b = X1[2] - X1[1]
	Y1_freq = (Y1 / float(a)) / float(b)
	os.mknod("histogram_X_axis.txt")
	np.savetxt("histogram_X_axis.txt", X1,fmt="%.5f",delimiter = '\t')
	os.mknod("histogram_Y_axis.txt")
	np.savetxt("histogram_Y_axis.txt", Y1_freq,fmt="%.5f",delimiter = '\t')
	# 直方图曲线
	Y3 = []
	#X3 = df_txt_r["value"]
	for i in range(len(X1)):
		y = (math.exp(-(X1[i]-mean) * (X1[i]-mean) / (2*std*std))) / (math.sqrt(2 * math.pi) * std)
		Y3.append(y)
	os.mknod("histogram_line_X.txt")
	np.savetxt("histogram_line_X.txt", X1,fmt="%.5f",delimiter = '\t')
	os.mknod("histogram_line_Y.txt")
	np.savetxt("histogram_line_Y.txt", Y3, fmt="%.5f", delimiter='\t')
	'''
	with open("histogram_line.txt",'w') as t:
		t.write(str(math.exp(-("x"-mean) * ("x"-mean) / (2*std*std))) / (math.sqrt(2 * math.pi) * std))
	'''
	
	# 输出QQ图数据
	yvals = np.arange(len(df_txt_r['value']))/float(len(df_txt_r['value']))
	Y2 = stats.norm.ppf(yvals)
	X2 = df_txt_r['value']
	os.mknod("QQ_Y_axis.txt")
	np.savetxt("QQ_Y_axis.txt", Y2[1:],fmt="%.5f",delimiter = '\t')
	os.mknod("QQ_X_axis.txt")
	f = open('QQ_X_axis.txt', 'w')
	X2[1:].to_csv('QQ_X_axis.txt', sep=' ', index=False)
	with open("QQ_line.txt",'w') as t:
		print >> t,'X1:%.5f\tY1:%.5f\nX2:%.5f\tY2:%.5f'% ((Y2[1]-abs(Y2[2]-Y2[1]))*std + mean,Y2[1]-abs(Y2[2]-Y2[1]),(Y2[n-1]+abs(Y2[n-1]-Y2[n-2]))*std + mean,Y2[n-1]+abs(Y2[n-1]-Y2[n-2]))
	data_output = open(out_file,'w')
	print >> data_output,'sample_num\taverage\tstdf\tskewness\tkurtosis\tks_D\tks_P\tsw_D\tsw_P'
	print >> data_output,'%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f' % (n,mean,std,skew,kurt,D1,P1,D2,P2)
	data_output.close()

def _main():
	parser = argparse.ArgumentParser(description='normal distance test')
	parser.add_argument('-i', '--input_data_file', help="input_data_file")
	parser.add_argument('-o', '--out_file', help="out_file")
	args = parser.parse_args()
	Normaldistest(args.input_data_file, args.out_file)

if __name__ == "__main__":
	_main()