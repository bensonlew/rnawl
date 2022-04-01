#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

import glob
import os


def make_estimators_table(input_path):
    os.chdir(input_path)
    file_path = glob.glob(r"{}/*.summary".format(input_path))
    # print(file_path)
    index_type = []
    sample_name = []
    line_list = []
    for fp in file_path:
        sp_name = fp.split('.')[-2]
        sample_name.append(sp_name)
        # print sample_name
        with open(fp, 'r') as f:
            index_type = f.readline().split('\t')
            index_type.pop(0)
            for line in f:
                line = line.strip().split('\t')
                line.pop(0)
                line = sp_name + '\t' + '\t'.join(line)
                line_list.append(line)
    out_file = os.path.join(input_path, 'estimators.xls')
    print(out_file)
    with open(out_file, 'w') as e:
        first_line = 'sample' + '\t' + '\t'.join(index_type)
        e.write(first_line)
        for line in line_list:
            e.write(line + '\n')
