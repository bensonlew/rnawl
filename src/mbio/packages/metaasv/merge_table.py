# -*- coding: utf-8 -*-

import sys
import logging
import os
import shutil
from Bio import SeqIO

input_dir = sys.argv[1]
output_path = sys.argv[2]
## table 获取二维表
out_asv_table = os.path.join(output_path, "ASV_table.xls")
listdirs = os.listdir(input_dir)
asv_md5_list = []
all_sample_abundance = {}
sample_list = []
for file in listdirs:
    file_path = os.path.join(input_dir, file)
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip().split("\t")
            if line[0] == "# Constructed from biom file":
                pass
            elif line[0] == "#OTU ID":
                sample = line[1]
                if sample not in sample_list:
                    sample_list.append(sample)
            elif line[0] == "ASV ID":
                sample = line[1]
                if sample not in sample_list:
                    sample_list.append(sample)
            else:
                asv_md5 = line[0]
                asc_abundance = line[1]
                if asv_md5 not in asv_md5_list:
                    asv_md5_list.append(asv_md5)
                if asv_md5 not in all_sample_abundance:
                    all_sample_abundance[asv_md5] = [{sample: asc_abundance}]
                else:
                    new_list = all_sample_abundance[asv_md5]
                    value = {sample: asc_abundance}
                    if value not in new_list:
                        new_list.append(value)
                    all_sample_abundance[asv_md5] = new_list


with open(out_asv_table, "w") as w:
    w.write("#OTU ID\t" + "\t".join(sample_list) + "\n")
    for md5_value in asv_md5_list:
        if md5_value in all_sample_abundance:
            w.write(md5_value + "\t")
            sample_abundance_list = all_sample_abundance[md5_value]
            n = 1
            for sample in sample_list:
                for new_sample_dict in sample_abundance_list:
                    if sample in new_sample_dict:
                        sample_abundance_value = new_sample_dict[sample]
                        if len(sample_list) == n:
                            w.write(sample_abundance_value)
                        else:
                            w.write(sample_abundance_value + "\t")
                if sample not in [x.keys()[0] for x in sample_abundance_list]:
                    print("sample: {}\tsample_abundance_list: {}".format(sample, [x.keys()[0] for x in sample_abundance_list]))
                    if len(sample_list) == n:
                        w.write("0.0")
                    else:
                        w.write("0.0" + "\t")
                n += 1
            w.write("\n")