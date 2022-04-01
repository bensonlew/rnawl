# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.24

import argparse


def slidling_window(s_step, input_file, out_file):
    """
    """
    chr_list = []
    win_data = {}
    with open(input_file, "r")as f, open(out_file, "w") as w:
        lines = f.readlines()
        for l in range(len(lines)):
            item = lines[l].strip().split("\t")
            if item[0] not in win_data.keys():
                chr_list.append(item[0])
                win_data[item[0]] = []
                i = 0
                x = 0
                step = s_step
            if int(item[1]) <= step:
                x += 1
            else:
                n = int(item[1]) / s_step
                win_data[item[0]].append(x)
                w.write(item[0] + "\t" + str(step) + "\t" + str(x) + "\n")
                i += 1
                if n > i:
                    m = i
                    for j in range(m, n):
                        step += s_step
                        x = 0
                        win_data[item[0]].append(x)
                        w.write(item[0] + "\t" + str(step) + "\t" + str(x) + "\n")
                        i += 1
                x = 1
                step += s_step
            if l < len(lines) - 1:
                tmp = lines[l+1].strip().split("\t")
                if item[0] != tmp[0]:
                    win_data[item[0]].append(x)
                    w.write(item[0] + "\t" + str(step) + "\t" + str(x) + "\n")
        try:
            win_data[item[0]].append(x)
            w.write(item[0] + "\t" + str(step) + "\t" + str(x) + "\n")
        except:
            print "{}文件为空！".format(input_file)


parser = argparse.ArgumentParser()
parser.add_argument("-step", "--step", help="step", required=True)
parser.add_argument("-i", "--input_file", help="input file", required=True)
parser.add_argument("-o", "--out_file", help="output file", required=True)
args = vars(parser.parse_args())


slidling_window(int(args["step"]), args["input_file"], args["out_file"])
