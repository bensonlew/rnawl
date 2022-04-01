# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# __version__ = 'v1.0'

import sys
import pandas as pd
import os
import random

def make_database(otu_table, outfile, index_num):
    """
    替换lefse分析的丰度表中的0的值,每一行的最小值
    """
    with open(otu_table, 'r') as f, open(outfile, 'w') as w:
        lines = f.readlines()
        w.write(lines[0])
        w.write(lines[1])
        if index_num == 3:
            w.write(lines[2])
        for line in lines[2:]:
            line = line.strip().split("\t")
            num_list = []
            if index_num in ['3', 3]:
                for num in line[3:]:
                    aa = float(num)
                    if aa != 0.0:
                       num_list.append(aa)
                if len(num_list) != 0:
                    min_num = min(num_list)
                    last_list = []
                    last_list.append(line[0])
                    for bb in line[1:]:
                        new_bb = float(bb)
                        if new_bb == 0.0:
                            new_cc = float(min_num / 10000)
                            last_list.append(new_cc)
                        else:
                            last_list.append(new_bb)
                    str_list = [str(i) for i in last_list]
                    w.write("\t".join(str_list) + "\n")
            elif index_num in ['2', 3]:
                for num in line[2:]:
                    aa = float(num)
                    if aa != 0.0:
                       num_list.append(aa)
                if len(num_list) != 0:
                    min_num = min(num_list)
                    last_list = []
                    last_list.append(line[0])
                    for bb in line[1:]:
                        new_bb = float(bb)
                        if new_bb == 0.0:
                            new_cc = float(min_num / 10000)
                            last_list.append(new_cc)
                        else:
                            last_list.append(new_bb)
                    str_list = [str(i) for i in last_list]
                    w.write("\t".join(str_list) + "\n")
            else:
                print("index_num不存在: %s"% index_num)


def replace_0_otu(otu_table,outfile, index_num):
    """
    lefse分析中的丰度表为0的数据用每一列的最小值的十万分之一进行替换
    """
    head = "head.xls"
    new_file = "out.xls"
    if index_num in ['3', 3]: ##因为
        os.system("head -n 2 {} > {}".format(otu_table, head))
        os.system("sed -e '1,2d' {} > {}".format(otu_table,new_file))
    elif index_num in ['2', 2]:
        os.system("head -n 1 {} > {}".format(otu_table, head))
        os.system("sed -e '1d' {} > {}".format(otu_table,new_file))
    else:
        os.system("head -n 1 {} > {}".format(otu_table, head))
        os.system("sed -e '1d' {} > {}".format(otu_table,new_file))
    data = pd.read_table(new_file, sep='\t', header=0)
    data_samples = []
    for i in data.columns:
        if i in ["#sample",'sample']:
            pass
        else:
            if i not in data_samples:
                data_samples.append(i)

    data = data.set_index("#sample") ##过滤掉样本总和都为0的结果
    data['Col_sum'] = data.apply(lambda x: x.sum(), axis=1)
    data = data[data['Col_sum'] != 0.0]
    del data['Col_sum']
    data = data.reset_index()

    for name in data_samples:##最小值进行替换
        sample_list = list(set(list(data[name].values)))
        # print(sample_list)
        sample_list.sort(reverse=True)
        sample_list.pop(-1)
        min_value = min(sample_list) / 10000
        random_value = float(random.uniform(0, min_value))/100
        new_value = random_value
        before_value = 0.0
        data[name].replace(before_value, new_value, inplace=True)
    data.to_csv(outfile, sep='\t', index=0)
    os.system("cat {} >> {}".format(outfile, head))
    os.system("mv {} {}".format(head, outfile))

if __name__ == "__main__":
    args = sys.argv[1:]
    otu_file = args[0]
    outfile = args[1]
    index_num = args[2]
    replace_0_otu(otu_file, outfile, index_num)
