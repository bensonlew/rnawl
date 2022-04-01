# # -*- coding: utf-8 -*-
# # __author__ = 'qindanhua'


def est_size(file_path):
    line_sum = []
    sample_num = 0
    with open(file_path, 'r') as f:
        f.readline()
        for lines in f:
            line = lines.strip().split('\t')
            line.pop(0)
            sample_num = len(line)
            ln_sum = sum(map(int, line))
            line_sum.append(ln_sum)
        size = sum(line_sum)/sample_num
        # print size
        return size

# est_size('in_otu_table.xls')
