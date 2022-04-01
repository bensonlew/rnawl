# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'


def reverse_table(input_table, output_table):
    write_line = []
    with open(input_table, "r") as f,open(output_table, "w") as w:
        sample = f.readline().strip().split("\t")[1:]
        write_line.append("sample")
        for s in sample:
            write_line.append(s)
        for line in f:
            line = line.strip().split("\t")
            for i in range(len(write_line)):
                write_line[i] += "\t{}".format(line[i])
        for wl in write_line:
            w.write("{}\n".format(wl))