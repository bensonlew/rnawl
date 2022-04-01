# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
import fileinput
import re
import os
import sys
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser(description="file")
parser.add_argument("-i", "--input", help="输入表格", required=True)
parser.add_argument("-method", "--method", help="row,column", required=True)
parser.add_argument("-combined", "--combined", help="合并小于此数值的区域", required=True)
parser.add_argument("-i1", "--middle_input", help="输入表格处理后的新输入文件，样本在行", required=True)
parser.add_argument("-i2", "--new_input", help="新表格转置后的文件，样本在列", required=True)
parser.add_argument("-o1", "--out1", help="输出的数值表", required=True)
parser.add_argument("-o2", "--out2", help="输出的百分表表格", required=True)
parser.add_argument("-group", "--group", help="分组", required=False)   # 第7个和第8个参数位置固定
parser.add_argument("-calculation ", "--calculation", help="sum,average,middle", required=False)
args = vars(parser.parse_args())

origin_input = args["input"]
method = args["method"]
combined = args["combined"]
middle_input = args["middle_input"]
new_input = args["new_input"]
out1 = args["out1"]
out2 = args["out2"]
if len(sys.argv) == 19:
    group = args["group"]
    calculation = args["calculation"]


def t_table(table_file, new_table):
    """
    转换颠倒表格内容
    """
    with open(table_file) as f, open(new_table, 'w') as w:
        lines = f.readlines()
        lines = [line for line in lines if (line != "\r\n") and (line != "\n")]
        lines = [line for line in lines if not re.search(r"^(\s*\t+?)\s*\t*\n*", line)]
        table_list = [i.rstrip().split('\t') for i in lines]
        table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
        w.writelines(table_list)


def group_table(input_file, group_table, output):
    """
    将根据分组文件，将输入文件进行组内合并，得到新的数据表格
    :param input_file: 输入的表格
    :param group_table: 分组文件
    :param output: 新的数据表格
    :return:
    """
    print "分组处理"
    sample_group = dict()  # 一个样本是属于哪个group的
    index_sample = dict()  # 一个OTU表中第几行属于哪个样本
    group_sample_num = defaultdict(int)  # 一个分组里面有多少的样本
    with open(group_table, "rb") as r:
        r.next()
        for line in r:
            line = line.rstrip().split("\t")
            sample_group[line[0]] = line[1]
            group_sample_num[line[1]] += 1

    with open(input_file, "rb") as r, open(output, 'wb') as w:
        group_list = list()
        for v in sample_group.values():
            group_list.append(v)
            group_list = list(set(group_list))

        line = r.next().rstrip().split("\t")
        for i in range(len(line)):
            index_sample[i] = line[i]

        w.write(index_sample[0] + "\t")
        w.write("\t".join(group_list) + "\n")
        for line in r:
            line = line.rstrip().split("\t")
            num = defaultdict(float)
            tmp = list()
            list1 = []
            mid_num = dict()

            w.write(line[0] + "\t")
            for i in range(1, len(line)):
                num[sample_group[index_sample[i]]] += round(float(line[i]),4)
            for m in group_list:
                for i in range(1, len(line)):
                    if sample_group[index_sample[i]] == m:
                        list1.append(float(line[i]))
                        if len(list1) == group_sample_num[m]:
                            list1.sort()
                            yu = float(group_sample_num[m]) % 2
                            index = int(float(group_sample_num[m]) / 2)
                            if yu == 0:
                                mid_num[m] = float(round(((float(list1[index - 1]) + float(list1[index])) / 2),4))
                                list1 = []
                            else:
                                mid_num[m] = list1[index]
                                list1 = []

            if calculation == "sum":
                for g in group_list:
                    tmp.append(str(num[g]))
            if calculation == "average":
                for g in group_list:
                    avg = float(round((num[g] / group_sample_num[g]),4))
                    tmp.append(str(avg))
            if calculation == "middle":
                for g in group_list:
                    tmp.append(str(mid_num[g]))
            w.write("\t".join(tmp))
            w.write("\n")


def get_final_result(input_file, combined, output, value):
    """
    将转置后的结果用来计算，得到柱形图画图数据和表格
    input_file:输入的表格
    combined：合并小于此数值的区域
    output：输出的百分比文件
    value：输出的丰度值文件
    :return:
    """
    with open(input_file) as fr, open(output, "w+")as fw, open(value, "w+")as fw2:
        lines = fr.readlines()
        lines = [line for line in lines if (line != "\r\n") and (line != "\n")]
        lines = [line for line in lines if not re.search(r"^(\s*\t+?)\s*\t*\n*", line)]
        first_line = lines[0].strip().split("\t")  # 物种名称列表
        # print first_line
        dic_origin = defaultdict(list)
        samples = []
        for line in lines[1:]:
            line_split = line.strip().split("\t")
            samples.append(line_split[0])  # 按顺序存放样本名称
            int_value = []
            for i in line_split[1:]:
                int_value.append(float(i))
            dic_origin[line_split[0]] = int_value
        # print dic_origin
        dic_percent = defaultdict(list)
        for key in dic_origin.keys():
            percent_list = []
            for percent_value in dic_origin[key]:
                if sum(dic_origin[key]) == 0:
                   raise Exception("{}对应的这组信息全为0，不能画图，请剔除，再做分析".format(key))
                else:
                    per = "%10f" % (percent_value / sum(dic_origin[key]))
                    percent_list.append(float(per))
            dic_percent[key] = percent_list
        # print dic_percent  # 百分比列表产生
        dic_compare = defaultdict(list)
        for i in range(len(first_line) - 1):
            for key in dic_percent.keys():
                dic_compare[i].append(dic_percent[key][i])
        # print dic_compare
        index = []
        for key in dic_compare.keys():
            put_or_not = []
            for per in dic_compare[key]:
                if float(per) > float(combined):
                    put_or_not.append(per)
                else:
                    pass
            if len(put_or_not) == 0:
                index.append(float(key))
        # print index
        new_names = []
        if len(index) != 0:  # 不为0，说明有需要合并的列
            new_dic = defaultdict(list)  # 处理之后的百分比表格
            for key in dic_percent.keys():
                value = []
                other = 0
                for i in range(len(first_line) - 1):
                    if i in index:
                        other += float(dic_percent[key][i])
                    else:
                        value.append(float(dic_percent[key][i]))
                value.append(other)
                new_dic[key] = value
            for i in range(len(first_line) - 1):
                if i not in index:
                    new_names.append(first_line[i + 1])
                else:
                    pass
            new_names.append("others")
            new_value_dic = defaultdict(list)  # 处理之后的数值表格
            for key in dic_origin.keys():
                value = []
                other = 0
                for i in range(len(first_line) - 1):
                    if i in index:
                        other += float(dic_origin[key][i])
                    else:
                        value.append(float(dic_origin[key][i]))
                value.append(other)
                new_value_dic[key] = value
        else:
            new_names = first_line[1:]
            new_dic = dic_percent
            new_value_dic = dic_origin
        new_first_line = "ID"
        for i in new_names:
            new_first_line = new_first_line + "\t" + str(i)
        fw.write(new_first_line + "\n")
        fw2.write(new_first_line + "\n")
        for j in samples:
            new_lines = j
            new_value_lines = j
            for value in new_dic[j]:
                new_lines += "\t" + str(value)
            fw.write(new_lines + "\n")
            for value2 in new_value_dic[j]:
                new_value_lines += "\t" + str(value2)
            fw2.write(new_value_lines + "\n")


def get_new_input_table(method, origin_input, middle_input, new_input):
    """
    通过行列转换，将输入的表格进行转置，统一样本名在第一行
    method:行列类型
    origin_input：原始的输入表格
    middle_input：样本在行的表格
    new_input：样本在列的表格
    :return:
    """
    if method == "row":  # 转换
        with open(origin_input, "r") as r, open(middle_input, "w+")as fw:
            lines = r.readlines()
            lines = [line for line in lines if (line != "\r\n") and (line != "\n")]
            lines = [line for line in lines if not re.search(r"^(\s*\t+?)\s*\t*\n*", line)]
            for line in lines:
                fw.write(line)
    if method == "column":  # 直接输出
        t_table(origin_input, middle_input)
    if len(sys.argv) == 19:  # 如果有分组，则组内合并
        group_table(middle_input, group, new_input)
    else:
        os.link(middle_input, new_input)
    newtable = new_input + '.T'
    t_table(new_input, newtable)
    get_final_result(newtable, combined, out1, out2)

get_new_input_table(method, origin_input, middle_input, new_input)




