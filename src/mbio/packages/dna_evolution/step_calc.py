# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20180907
from sys import argv
import argparse
import os
import json
parse = argparse.ArgumentParser(description="用于变异位点比较分析染色体分布图滑窗大小传入")
parse.add_argument("-p", "--path", help="传入文件路径，如pop.table，或者win.stat", required=True)
parse.add_argument("-s", "--step_num", help="输入步长", required=True)
parse.add_argument("-o", "--output", help="输出结果", required=True)
parse.add_argument("-v", "--variant_type", help="变异类型", required=True)
args = vars(parse.parse_args())
result_path = args["path"]
output_path = args["output"]
variant_type = args["variant_type"]
s_step = int(args["step_num"])


def step_calc(variants_type):
    chr_lists = []
    win_data = {}
    step = ''
    x = 0
    i = 0
    # 最后与图对应的信息为win_data[chr1] = [5,8,9,2,7,12]
    # x 表示的是list中的元素，i记录的为相当于list的index的值，下面当n > i 时的循环用于将步长内为空的标记为0。
    # tmp 用于判断将不同的染色体归于同一类。
    if variants_type == "snp":
        with open(os.path.join(output_path, "file_snp.txt"), "r")as f:
            lines = f.readlines()
            chr_lists = []  # 用于存放染色体
            win_data = {}  # 用于存放染色体内分段的个数
            i = 0  # 定义存放在list中第一个位置的数值个数
            step = 0  # 在全局可能使用，所以定义一个全局变量
            for l in range(len(lines)):
                item = lines[l].strip().split("\t")
                if item[0] not in win_data:
                    win_data[item[0]] = []
                    chr_lists.append(item[0])
                if int(item[1]) >= int(step):
                    if i != 0:
                        win_data[item[0]].append(i)
                        step = step + s_step
                        i = 0
                        if int(item[1]) < int(step):
                            i = i + 1
                        else:
                            m = (int(item[1]) - int(step)) / int(s_step) + 1
                            for m in range(m):
                                win_data[item[0]].append(0)
                                step = step + s_step
                            if int(item[1]) < int(step):
                                i = i + 1
                    else:
                        n = (int(item[1]) - int(step)) / int(s_step)
                        step = s_step
                        for j in range(n):
                            win_data[item[0]].append(0)
                            step = step + s_step
                        if int(item[1]) < int(step):
                            i = i + 1
                else:
                    i = i + 1
                if l < len(lines) - 1:
                    temp = lines[l + 1].strip().split("\t")[0]
                    if item[0] != temp:
                        win_data[item[0]].append(i)
                        i = 0
                        step = 0
                else:
                    win_data[item[0]].append(i)
        s_max, start, end = 0, 0, 0
        chr_list, sca_list, other_list, final_list = [], [], [], []
        for i in chr_lists:  # 对染色体排序
            if i.startswith("chr"):
                num = i.strip().split("chr")[-1]
                chr_list.append(int(num))
            elif i.startswith("sca"):
                num = i.strip().split("sca")[-1]
                sca_list.append(int(num))
            else:
                other_list.append(i)
        chr_list.sort()
        sca_list.sort()
        for i in chr_list:
            final_list.append("chr" + str(i))
        for i in sca_list:
            final_list.append("sca" + str(i))
        for i in other_list:
            final_list.append(i)
        for chr in final_list:
            s_max = len(win_data[chr]) * s_step
            if s_max > end:
                end = s_max
        if os.path.exists(os.path.join(output_path, "data_snp")):
            os.remove(os.path.join(output_path, "data_snp"))
        with open (os.path.join(output_path, "data_snp"), "w") as f:
            f.write(",".join(final_list) + "\n")
            f.write(json.dumps(win_data) + "\n")
            f.write(str(start) + "\n")
            f.write(str(end))
    elif variants_type == "indel":
        with open(os.path.join(output_path, "file_indel.txt"), "r")as f:
            lines = f.readlines()
            chr_lists = []  # 用于存放染色体
            win_data = {}  # 用于存放染色体内分段的个数
            i = 0  # 定义存放在list中第一个位置的数值个数
            step = 0  # 在全局可能使用，所以定义一个全局变量
            for l in range(len(lines)):
                item = lines[l].strip().split("\t")
                if item[0] not in win_data:
                    win_data[item[0]] = []
                    chr_lists.append(item[0])
                if int(item[1]) >= int(step):
                    if i != 0:
                        win_data[item[0]].append(i)
                        step = step + s_step
                        i = 0
                        if int(item[1]) < int(step):
                            i = i + 1
                        else:
                            m = (int(item[1]) - int(step)) / int(s_step) + 1
                            for m in range(m):
                                win_data[item[0]].append(0)
                                step = step + s_step
                            if int(item[1]) < int(step):
                                i = i + 1
                    else:
                        n = (int(item[1]) - int(step)) / int(s_step)
                        step = s_step
                        for j in range(n):
                            win_data[item[0]].append(0)
                            step = step + s_step
                        if int(item[1]) < int(step):
                            i = i + 1
                else:
                    i = i + 1
                if l < len(lines) - 1:
                    temp = lines[l + 1].strip().split("\t")[0]
                    if item[0] != temp:
                        win_data[item[0]].append(i)
                        i = 0
                        step = 0
                else:
                    win_data[item[0]].append(i)
        s_max, start, end = 0, 0, 0
        chr_list, sca_list, other_list, final_list = [], [], [], []
        for i in chr_lists:  # 对染色体排序
            if i.startswith("chr"):
                num = i.strip().split("chr")[-1]
                chr_list.append(int(num))
            elif i.startswith("sca"):
                num = i.strip().split("sca")[-1]
                sca_list.append(int(num))
            else:
                other_list.append(i)
        chr_list.sort()
        sca_list.sort()
        for i in chr_list:
            final_list.append("chr" + str(i))
        for i in sca_list:
            final_list.append("sca" + str(i))
        for i in other_list:
            final_list.append(i)
        for chr in final_list:
            s_max = len(win_data[chr]) * s_step
            if s_max > end:
                end = s_max
        if os.path.exists(os.path.join(output_path, "data_indel")):
            os.remove(os.path.join(output_path, "data_indel"))
        with open(os.path.join(output_path, "data_indel"), "w") as f:
            f.write(",".join(final_list) + "\n")
            f.write(json.dumps(win_data) + "\n")
            f.write(str(start) + "\n")
            f.write(str(end))
    elif variants_type == "all":
        with open(os.path.join(output_path, "file_all.txt"), "r")as f:
            lines = f.readlines()
            chr_lists = []  # 用于存放染色体
            win_data = {}  # 用于存放染色体内分段的个数
            i = 0  # 定义存放在list中第一个位置的数值个数
            step = 0  # 在全局可能使用，所以定义一个全局变量
            for l in range(len(lines)):
                item = lines[l].strip().split("\t")
                if item[0] not in win_data:
                    win_data[item[0]] = []
                    chr_lists.append(item[0])
                if int(item[1]) >= int(step):
                    if i != 0:
                        win_data[item[0]].append(i)
                        step = step + s_step
                        i = 0
                        if int(item[1]) < int(step):
                            i = i + 1
                        else:
                            m = (int(item[1]) - int(step)) / int(s_step) + 1
                            for m in range(m):
                                win_data[item[0]].append(0)
                                step = step + s_step
                            if int(item[1]) < int(step):
                                i = i + 1
                    else:
                        n = (int(item[1]) - int(step)) / int(s_step)
                        step = s_step
                        for j in range(n):
                            win_data[item[0]].append(0)
                            step = step + s_step
                        if int(item[1]) < int(step):
                            i = i + 1
                else:
                    i = i + 1
                if l < len(lines) - 1:
                    temp = lines[l + 1].strip().split("\t")[0]
                    if item[0] != temp:
                        win_data[item[0]].append(i)
                        i = 0
                        step = 0
                else:
                    win_data[item[0]].append(i)
        # with open(os.path.join(output_path, "file_all.txt"), "r")as f:
        #     lines = f.readlines()
        #     for l in range(len(lines)):
        #         item = lines[l].strip().split("\t")
        #         if item[0] not in win_data.keys():
        #             win_data[item[0]] = []
        #             chr_lists.append(item[0])
        #             i = 0
        #             x = 0
        #             step = s_step
        #         if int(item[1]) <= step:
        #             x += 1
        #         else:
        #             n = int(item[1]) / s_step
        #             win_data[item[0]].append(x)
        #             i += 1
        #             if n > i:
        #                 m = i
        #                 for j in range(m, n):
        #                     step += s_step
        #                     x = 0
        #                     win_data[item[0]].append(x)
        #                     i += 1
        #             x = 1
        #             step += s_step
        #         if l < len(lines) - 1:
        #             tmp = lines[l + 1].strip().split("\t")
        #             if item[0] != tmp[0]:
        #                 win_data[item[0]].append(x)
        s_max, start, end = 0, 0, 0
        chr_list, sca_list, other_list, final_list = [], [], [], []
        for i in chr_lists:  # 对染色体排序
            if i.startswith("chr"):
                num = i.strip().split("chr")[-1]
                chr_list.append(int(num))
            elif i.startswith("sca"):
                num = i.strip().split("sca")[-1]
                sca_list.append(int(num))
            else:
                other_list.append(i)
        chr_list.sort()
        sca_list.sort()
        for i in chr_list:
            final_list.append("chr" + str(i))
        for i in sca_list:
            final_list.append("sca" + str(i))
        for i in other_list:
            final_list.append(i)
        for chr in final_list:
            s_max = len(win_data[chr]) * s_step
            if s_max > end:
                end = s_max
        if os.path.exists(os.path.join(output_path, "data_all")):
            os.remove(os.path.join(output_path, "data_all"))
        with open(os.path.join(output_path, "data_all"), "w") as f:
            f.write(",".join(final_list) + "\n")
            f.write(json.dumps(win_data) + "\n")
            f.write(str(start) + "\n")
            f.write(str(end))


def get_file():
    chr_list_snp = []
    pos_list_snp = []
    chr_list_indel = []
    pos_list_indel = []
    chr_list_all = []
    pos_list_all = []
    with open(result_path, "r") as r:
        lines = r.readlines()
        for line in lines[1:]:
            rows = line.strip().split()
            chr_list_all.append(rows[0])
            pos_list_all.append(rows[1])
            if rows[2] == "SNP":
                chr_list_snp.append(rows[0])
                pos_list_snp.append(rows[1])
            elif rows[2] == "INDEL":
                chr_list_indel.append(rows[0])
                pos_list_indel.append(rows[1])
            else:
                raise Exception("基因型格式不正确")

    with open(os.path.join(args["output"], "file_snp.txt"), "w") as w:
        for (i, j) in zip(chr_list_snp, pos_list_snp):
            w.write(i + "\t" + j + "\n")
    with open(os.path.join(args["output"], "file_indel.txt"), "w") as w:
        for (i, j) in zip(chr_list_indel, pos_list_indel):
            w.write(i + "\t" + j + "\n")
    with open(os.path.join(args["output"], "file_all.txt"), "w") as w:
        for (i, j) in zip(chr_list_all, pos_list_all):
            w.write(i + "\t" + j + "\n")

get_file()
step_calc(variant_type)


