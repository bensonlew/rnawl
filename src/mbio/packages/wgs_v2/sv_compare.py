# -*- coding: utf-8 -*-
# __author__: zengjing

import re
import os
import json
import argparse


def region_filter(item, region_info):
    """
    区域筛选
    """
    if item[0] not in region_info.keys():
        return False
    for start_end in region_info[item[0]]:
        start_end = start_end.split("-")
        if start_end[0] and start_end[1]:
            if int(start_end[0]) <= int(item[1]) <= int(start_end[1]):
                return True
        elif start_end[0]:
            if int(item[1]) >= int(start_end[0]):
                return True
        elif start_end[1]:
            if int(item[1]) <= int(start_end[1]):
                return True
        else:
            return True
    return False

def single_sv_compare(sv_file, config_file, output_dir):
    """
    单条件样本间SV比较
    sv_file：样本的sv vcf文件
    config_file：筛选条件配置文件,json格式
    output_dir: 输出结果文件夹
    """
    sample_info = json.loads(open(config_file).read())
    if "Region" in sample_info.keys():
        region_info = sample_info["Region"]
        sample_info.pop("Region")
    sample_list, file_info, total_info, cs_chr = [], {}, {}, {}
    for cs in sample_info.keys():
        sample = cs.split("_VS_")[0]
        c_sample = cs.split("_VS_")[1]
        total_info[cs], cs_chr[cs] = {}, []
        sv_detail = os.path.join(output_dir, cs + ".detail.xls")
        sv_total = os.path.join(output_dir, cs + ".summary.xls")
        file_info[cs] = {"detail": open(sv_detail, "w"), "total": open(sv_total, "w")}
        file_info[cs]["detail"].write("#SV ID\tCHROM \tStart\tEnd\tLength\t{}_genotype\t{}_depth\t{}_genotype\t{}_depth\n".format(sample, sample, c_sample, c_sample))
        file_info[cs]["total"].write("#CHROM\tDEL\tINS\tDUP\tINV\tBND\n")
    sample_index, type_index = {}, {}
    with open(sv_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("##"):
                continue
            item = line.strip().split("\t")
            if line.startswith("#CHROM"):  # 样本的index
                for i in range(9, len(item)):
                    sample_index[item[i]] = i
            if item[6] != "PASS":  # 过滤模式不为PASS则该位点是SV的概率低，因此过滤掉
                continue
            if region_info != "all":  # 区域过滤
                status = region_filter(item, region_info)
                if not status:
                    continue
            ty = item[8].split(":")  # 确定每种类型所在的位置
            for i in range(len(ty)):
                type_index[ty[i]] = i
            m = re.match(r".+;SVTYPE=(.+?);.+END=(\d+);.+", item[7])  # 计算长度
            sv_type = m.group(1)
            end = m.group(2)
            if sv_type in ["DEL", "DUP", "INV"]:
                length = int(end) - int(item[1])
            elif sv_type == "INS":
                n = re.match(r".+;INSLEN=(\d+);.+", item[7])
                length = n.group(1)
            elif sv_type == "BND":
                length = len(item[3])
            else:
                continue
            for cs in sample_info.keys():  # 每个比较样本比较组进行过滤
                sample = cs.split("_VS_")[0]
                compare_sample = cs.split("_VS_")[1]
                sample_item = item[sample_index[sample]].split(":")
                compare_sample_item = item[sample_index[compare_sample]].split(":")
                compare_type = sample_info[cs]["compare_type"]  # 过滤两样本间相同/不同/全部的位点
                s_ft = sample_item[type_index["FT"]]
                cs_ft = compare_sample_item[type_index["FT"]]
                if s_ft != "PASS" and cs_ft != "PASS":
                    continue
                if compare_type == "same":
                    if s_ft != "PASS" or cs_ft != "PASS":
                        continue
                elif compare_type == "diff":
                    if s_ft == "PASS" and cs_ft == "PASS":
                        continue
                else:
                    pass
                depth = int(sample_item[type_index["DR"]]) + int(sample_item[type_index["DV"]])
                c_depth = int(compare_sample_item[type_index["DR"]]) + int(compare_sample_item[type_index["DV"]])
                if sample in sample_info[cs].keys():
                    if "geno_type" in sample_info[cs][sample].keys():
                        gt = list(set(sample_item[type_index["GT"]].split("/")))
                        if sample_info[cs][sample]["geno_type"] == "same":
                            if len(gt) != 1:
                                continue
                        elif sample_info[cs][sample]["geno_type"] == "mixed":
                            if len(gt) != 2:
                                continue
                    if "depth" in sample_info[cs][sample].keys():
                        min = sample_info[cs][sample]["depth"].split("-")[0]
                        max = sample_info[cs][sample]["depth"].split("-")[1]
                        if min:
                            if depth < int(min):
                                continue
                        if max:
                            if depth > int(max):
                                continue
                if compare_sample in sample_info[cs].keys():
                    if "geno_type" in sample_info[cs][compare_sample].keys():
                        gt = list(set(compare_sample_item[type_index["GT"]].split("/")))
                        if sample_info[cs][compare_sample]["geno_type"] == "same":
                            if len(gt) != 1:
                                continue
                        elif sample_info[cs][compare_sample]["geno_type"] == "mixed":
                            if len(gt) != 2:
                                continue
                    if "depth" in sample_info[cs][compare_sample].keys():
                        min = sample_info[cs][compare_sample]["depth"].split("-")[0]
                        max = sample_info[cs][compare_sample]["depth"].split("-")[1]
                        if min:
                            if c_depth < int(min):
                                continue
                        if max:
                            if c_depth > int(max):
                                continue
                sample_genotype = sv_type if s_ft == "PASS" else "-"
                c_sample_genotype = sv_type if cs_ft == "PASS" else "-"
                if item[0] not in cs_chr[cs]:
                    cs_chr[cs].append(item[0])
                    total_info[cs][item[0]] = {"DEL": 0, "INS": 0, "DUP": 0, "INV": 0, "BND": 0}
                total_info[cs][item[0]][sv_type] += 1
                file_info[cs]["detail"].write(item[2] + "\t" + item[0] + "\t" + item[1] + "\t" + end + "\t" + str(length) + "\t")
                file_info[cs]["detail"].write(sample_genotype + "\t" + str(depth) + "\t" + c_sample_genotype + "\t" + str(c_depth) + "\n")
    for cs in sample_info.keys():
        file_info[cs]["detail"].close()
        for chr in cs_chr[cs]:
            file_info[cs]["total"].write(chr + "\t" + str(total_info[cs][chr]["DEL"]) + "\t" + str(total_info[cs][chr]["INS"]) + "\t")
            file_info[cs]["total"].write(str(total_info[cs][chr]["DUP"]) + "\t" + str(total_info[cs][chr]["INV"]) + "\t")
            file_info[cs]["total"].write(str(total_info[cs][chr]["BND"]) + "\n")
        file_info[cs]["total"].close()

def multiple_sv_compare(sv_file, config_file, output_dir):
    """
    多条件样本间SV比较
    sv_file：样本的sv vcf文件
    config_file：筛选条件配置文件,json格式
    output_dir: 输出结果文件夹
    """
    sample_info = json.loads(open(config_file).read())
    if "Region" in sample_info.keys():
        region_info = sample_info["Region"]
        sample_info.pop("Region")
    sample_list, chr_list_, file_info, total_info = [], [], {}, {}
    sv_detail = open(os.path.join(output_dir, "sv.detail.xls"), "w")
    sv_total = open(os.path.join(output_dir, "sv.summary.xls"), "w")
    sv_detail.write("#SV ID\tCHROM \tStart\tEnd\tLength")
    sv_total.write("#CHROM\tDEL\tINS\tDUP\tINV\tBND\n")
    for cs in sample_info.keys():
        sample = cs.split("_VS_")[0]
        c_sample = cs.split("_VS_")[1]
        if sample not in sample_list:
            sample_list.append(sample)
        if c_sample not in sample_list:
            sample_list.append(c_sample)
    for sample in sample_list:
        sv_detail.write("\t{}_genotype\t{}_depth".format(sample, sample))
    sv_detail.write("\n")
    sample_index, type_index = {}, {}
    with open(sv_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("##"):
                continue
            item = line.strip().split("\t")
            if line.startswith("#CHROM"):  # 样本的index
                for i in range(9, len(item)):
                    sample_index[item[i]] = i
            if item[6] != "PASS":  # 过滤模式不为PASS则该位点是SV的概率低，因此过滤掉
                continue
            if region_info != "all":  # 区域过滤
                status = region_filter(item, region_info)
                if not status:
                    continue
            ty = item[8].split(":")  # 确定每种类型所在的位置
            for i in range(len(ty)):
                type_index[ty[i]] = i
            ft_list, break_flag = [], False
            for cs in sample_info.keys():  # 每个比较样本比较组进行过滤
                sample = cs.split("_VS_")[0]
                compare_sample = cs.split("_VS_")[1]
                sample_item = item[sample_index[sample]].split(":")
                compare_sample_item = item[sample_index[compare_sample]].split(":")
                compare_type = sample_info[cs]["compare_type"]  # 过滤两样本间相同/不同/全部的位点
                s_ft = sample_item[type_index["FT"]]
                cs_ft = compare_sample_item[type_index["FT"]]
                ft_list.append(s_ft)
                ft_list.append(cs_ft)
                if s_ft != "PASS" and cs_ft != "PASS":
                    break_flag = True
                    break
                if compare_type == "same":
                    if s_ft != "PASS" or cs_ft != "PASS":
                        break_flag = True
                        break
                elif compare_type == "diff":
                    if s_ft == "PASS" and cs_ft == "PASS":
                        break_flag = True
                        break
                else:
                    pass
                depth = int(sample_item[type_index["DR"]]) + int(sample_item[type_index["DV"]])
                c_depth = int(compare_sample_item[type_index["DR"]]) + int(compare_sample_item[type_index["DV"]])
                if sample in sample_info[cs].keys():
                    if "geno_type" in sample_info[cs][sample].keys():
                        gt = list(set(sample_item[type_index["GT"]].split("/")))
                        if sample_info[cs][sample]["geno_type"] == "same":
                            if len(gt) != 1:
                                break_flag = True
                                break
                        elif sample_info[cs][sample]["geno_type"] == "mixed":
                            if len(gt) != 2:
                                break_flag = True
                                break
                    if "depth" in sample_info[cs][sample].keys():
                        min = sample_info[cs][sample]["depth"].split("-")[0]
                        max = sample_info[cs][sample]["depth"].split("-")[1]
                        if min:
                            if depth < int(min):
                                break_flag = True
                                break
                        if max:
                            if depth > int(max):
                                break_flag = True
                                break
                if compare_sample in sample_info[cs].keys():
                    if "geno_type" in sample_info[cs][compare_sample].keys():
                        gt = list(set(compare_sample_item[type_index["GT"]].split("/")))
                        if sample_info[cs][compare_sample]["geno_type"] == "same":
                            if len(gt) != 1:
                                break_flag = True
                                break
                        elif sample_info[cs][compare_sample]["geno_type"] == "mixed":
                            if len(gt) != 2:
                                break_flag = True
                                break
                    if "depth" in sample_info[cs][compare_sample].keys():
                        min = sample_info[cs][compare_sample]["depth"].split("-")[0]
                        max = sample_info[cs][compare_sample]["depth"].split("-")[1]
                        if min:
                            if c_depth < int(min):
                                break_flag = True
                                break
                        if max:
                            if c_depth > int(max):
                                break_flag = True
                                break
            if break_flag:
                continue
            if "PASS" not in ft_list:  # 如果所有样本的FT都没有PASS则在这次比较中将这个SV过滤掉
                continue
            # print item[:3]
            m = re.match(r".+;SVTYPE=(.+?);.+END=(\d+);.+", item[7])  # 计算长度
            sv_type = m.group(1)
            end = m.group(2)
            sample_genotype = sv_type if s_ft == "PASS" else "-"
            c_sample_genotype = sv_type if cs_ft == "PASS" else "-"
            if sv_type in ["DEL", "DUP", "INV"]:
                length = int(end) - int(item[1])
            elif sv_type == "INS":
                n = re.match(r".+;INSLEN=(\d+);.+", item[7])
                if not n:
                    continue
                length = n.group(1)
                end = int(item[1]) + int(length)
            elif sv_type == "BND":
                length = len(item[3])
                end = int(item[1]) + length
            else:
                continue
            if item[0] not in chr_list_:
                chr_list_.append(item[0])
                total_info[item[0]] = {"DEL": 0, "INS": 0, "DUP": 0, "INV": 0, "BND": 0}
            total_info[item[0]][sv_type] += 1
            sv_detail.write(item[2] + "\t" + item[0] + "\t" + item[1] + "\t" + end + "\t" + str(length))
            for sample in sample_list:
                sample_item = item[sample_index[sample]].split(":")
                ft = sample_item[type_index["FT"]]
                sample_genotype = sv_type if ft == "PASS" else "-"
                depth = int(sample_item[type_index["DR"]]) + int(sample_item[type_index["DV"]])
                sv_detail.write("\t" + sample_genotype + "\t" + str(depth))
            sv_detail.write("\n")
    sv_detail.close()
    for cs in sample_info.keys():
        for chr in chr_list_:
            sv_total.write(chr + "\t" + str(total_info[chr]["DEL"]) + "\t" + str(total_info[chr]["INS"]) + "\t")
            sv_total.write(str(total_info[chr]["DUP"]) + "\t" + str(total_info[chr]["INV"]) + "\t")
            sv_total.write(str(total_info[chr]["BND"]) + "\n")
    sv_total.close()


# params = {'Region': {'chr2': '2-1000000', 'chr1': '1-'}, 'SRR5739119_VS_SRR5739120': {'SRR5739120': {'geno_type': 'mixed', 'depth': '0-20'}, 'compare_type': 'same', 'SRR5739119': {'geno_type': 'same', 'depth': '5-10'}}}

parser = argparse.ArgumentParser(description='SV两两样本间进行比较')
parser.add_argument('-sv', '--sv_file', help='file:sv.vcf', required=True)
parser.add_argument('-c', '--config_file', help='params of config file', required=True)
parser.add_argument('-o', '--output_dir', help='output dir', required=True)
parser.add_argument('-m', '--model', help='single or multiple', required=True)

args = vars(parser.parse_args())
if args["model"] == "single":
    single_sv_compare(args['sv_file'], args['config_file'], args['output_dir'])
elif args["model"] == "multiple":
    multiple_sv_compare(args['sv_file'], args['config_file'], args['output_dir'])
else:
    raise Exception("model模式只能是single或multiple")
