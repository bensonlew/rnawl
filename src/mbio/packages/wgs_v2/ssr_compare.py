# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190328

import os
import re
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

def get_repeat_count(item):
    """
    求样本对应的等位基因单元重复数
    """
    gene = item.split(":")[0]
    count = item.split(":")[1]
    sum = 0
    for i in count.split(","):
        sum += int(i) * len(gene)
    return sum

def diff_double_list(a_list, b_list):
    """
    求两个list之间的差集
    """
    coll = []
    for i in a_list:
        if i not in b_list:
            coll.append(i)
    for i in b_list:
        if i not in a_list:
            coll.append(i)
    return coll

def sample_filter(sample_info, sample_index, item):
    """
    样本基因型、深度、比较类型过滤,返回过滤状态,每个样本的重复单元数和测序深度
    """
    status = True
    sample_stat = {}
    for cs in sample_info.keys():
        c_sample = cs.split("_VS_")[0]
        s_sample = cs.split("_VS_")[1]
        c_gene_type = list(set(item[sample_index[c_sample]].split(":")[1].split(",")))
        s_gene_type = list(set(item[sample_index[s_sample]].split(":")[1].split(",")))
        coll = diff_double_list(c_gene_type, s_gene_type)
        if sample_info[cs][c_sample]["gene_type"] == "same" and len(c_gene_type) != 1:
            status = False
            break
        elif sample_info[cs][c_sample]["gene_type"] == "mixed" and len(c_gene_type) == 1:
            status = False
            break
        if sample_info[cs][s_sample]["gene_type"] == "same" and len(s_gene_type) != 1:
            status = False
            break
        elif sample_info[cs][s_sample]["gene_type"] == "mixed" and len(s_gene_type) == 1:
            status = False
            break
        c_count = get_repeat_count(item[sample_index[c_sample]]) if c_gene_type != ["-"] else "-"
        s_count = get_repeat_count(item[sample_index[s_sample]]) if s_gene_type != ["-"] else "-"
        compare_type = sample_info[cs]["compare_type"]  # 样本间是否相同过滤
        if compare_type == "same" and coll:
            status = False
            break
        if compare_type == "diff" and not coll:
            status = False
            break
        c_depth = item[sample_index[c_sample]].split(":")[-1]  # 测序深度过滤
        s_depth = item[sample_index[s_sample]].split(":")[-1]
        c_min_depth = sample_info[cs][c_sample]["depth"].split("-")[0]
        c_max_depth = sample_info[cs][c_sample]["depth"].split("-")[1]
        s_min_depth = sample_info[cs][s_sample]["depth"].split("-")[0]
        s_max_depth = sample_info[cs][s_sample]["depth"].split("-")[1]
        if c_depth != "-":
            if c_min_depth and int(c_depth) < int(c_min_depth):
                status = False
                break
            if c_max_depth and int(c_depth) > int(c_max_depth):
                status = False
                break
        if s_depth != "-":
            if s_min_depth and int(s_depth) < int(s_min_depth):
                status = False
                break
            if s_max_depth and int(s_depth) > int(s_max_depth):
                status = False
                break
        sample_stat[c_sample] = [c_count, c_depth]
        sample_stat[s_sample] = [s_count, s_depth]
    if status:
        return True, sample_stat
    else:
        return False, sample_stat


def group_filter(group_info, sample_index, item):
    """
    组的深度、基因型频率、缺失率过滤，返回过滤状态，基因型频率和缺失率值
    """
    group_stat = {}  #
    status, group_miss_fre, group_gene_fre = True, 0, 0
    for group in group_info.keys():
        g_min_depth = group_info[group]["depth"].split("-")[0]
        g_max_depth = group_info[group]["depth"].split("-")[1]
        g_min_gene_fre = group_info[group]["gene_fre"].split("-")[0]
        g_max_gene_fre = group_info[group]["gene_fre"].split("-")[1]
        g_miss_fre = float(group_info[group]["miss_fre"])
        group_sample = group_info[group]["sample"]
        miss_value, gene_value = 0, 0
        for sample in group_sample:
            s_depth = item[sample_index[sample]].split(":")[-1]
            s_gene_type = list(set(item[sample_index[sample]].split(":")[1].split(",")))
            s_count = get_repeat_count(item[sample_index[sample]]) if s_gene_type != ["-"] else "-"  # 样本的重复单元长度总和
            if s_count == "-":  # 没有这个SSR此处算缺失
                miss_value += 1
            elif g_min_depth and int(s_depth) < int(g_min_depth):  # 不符合深度的算缺失
                miss_value += 1
            elif g_max_depth and int(s_depth) > int(g_max_depth):
                miss_value += 1
            elif s_count == int(item[5]):  # 样本的重复单元长度和参考基因组一致时，算基因型频率
                gene_value += 1
        group_miss_fre = float(miss_value) / len(group_sample)  # 缺失率过滤
        group_gene_fre = float(gene_value) / len(group_sample)  # 基因型频率过滤
        if g_miss_fre < group_miss_fre:
            status = False
            break
        if g_min_gene_fre and float(g_min_gene_fre) > group_gene_fre:
            status = False
            break
        elif g_max_gene_fre and float(g_max_gene_fre) < group_gene_fre:
            status = False
            break
        group_stat[group] = [group_miss_fre, group_gene_fre]
    return status, group_stat


def single_ssr_compare(ssr_file, config_file, output_dir):
    """
    SSR单条件样本间批量比较
    params: {"region_select": {"chr1": 1-10000}, "allele_num": 2/3, "compare_type": "sample/diff",
            "gene_type": "same/mixed", "depth": "2-100", "sample_info": [sample1_VS_sample2)]}
            allele_num: 等位基因数目(2/3)
            compare_type: 比较类型(same/diff)
            gene_type: 基因型(same/mixed)
            depth: 测序深度(-分隔)
            sample_info: 进行比较的样本([(sample1,sample2),(sample2,sample3)])
    """
    params = json.loads(open(config_file).read())
    region_select = params["region_select"]
    compare_type = params["compare_type"]
    gene_type = params["gene_type"]
    sample_info = params["sample_info"]
    min_depth = params["depth"].split("-")[0]
    max_depth = params["depth"].split("-")[1]
    file_info, stat_info, sample_index = {}, {}, {}
    for cs in sample_info:
        stat_info[cs] = {}
        ssr_detail = os.path.join(output_dir, cs + ".detail.xls")
        ssr_stat = os.path.join(output_dir, cs + ".stat.xls")
        file_info[cs] = {"detail": open(ssr_detail, "w"), "stat": open(ssr_stat, "w")}
        file_info[cs]["detail"].write("#CHROM\tStart\tEnd\tSSR ID\tRepeat unit\tReference Repeat count\t")
        file_info[cs]["detail"].write("{}_repeat_count\t{}_allele_depth\t".format(cs.split("_VS_")[0], cs.split("_VS_")[0]))
        file_info[cs]["detail"].write("{}_repeat_count\t{}_allele_depth\n".format(cs.split("_VS_")[1], cs.split("_VS_")[1]))
        file_info[cs]["stat"].write("#CHROM\tSSR Number\tc\tc*\tp1\tp2\tp3\tp4\tp5\tp6\n")
    with open(ssr_file, "rb") as f:
        for line in f:
            if line.startswith("##"):
                continue
            item = line.strip().split("\t")
            if line.startswith("#CHROM"):  # 样本的index
                for i in range(9, len(item)):
                    sample_index[item[i]] = i
                continue
            # if item[6] != "PASS":
            #     continue
            if region_select != "all":  # 区域过滤
                status = region_filter(item, region_select)
                if not status:
                    continue
            # print item
            # ref_gene_type = list(set(re.findall("\d+", item[4])))
            # if len(ref_gene_type) not in allele_num:  # 等位基因数目过滤(过滤参考基因组)
            #     continue
            for cs in sample_info:  # 纯合杂合过滤
                c_sample = cs.split("_VS_")[0]
                s_sample = cs.split("_VS_")[1]
                print sample_index
                c_gene_type = list(set(item[sample_index[c_sample]].split(":")[1].split(",")))
                s_gene_type = list(set(item[sample_index[s_sample]].split(":")[1].split(",")))
                if c_gene_type == ["-"] and s_gene_type == ["-"]:
                    continue
                diff_coll = diff_double_list(c_gene_type, s_gene_type)  # 相同不相同过滤
                if diff_coll and compare_type == "same":
                    continue
                if not diff_coll and compare_type == "diff":
                    continue
                if gene_type == "same":
                    if len(c_gene_type) != 1 or len(s_gene_type) != 1:
                        continue
                elif gene_type == "mixed":
                    if len(c_gene_type) == 1 or len(s_gene_type) == 1:
                        continue
                c_count = get_repeat_count(item[sample_index[c_sample]]) if c_gene_type != ["-"] else "-"
                s_count = get_repeat_count(item[sample_index[s_sample]]) if s_gene_type != ["-"] else "-"
                c_depth = item[sample_index[c_sample]].split(":")[-1]  # 测序深度过滤
                s_depth = item[sample_index[s_sample]].split(":")[-1]
                if c_depth != "-":
                    if min_depth and int(c_depth) < int(min_depth):
                        continue
                    if max_depth and int(c_depth) > int(max_depth):
                        continue
                if s_depth != "-":
                    if min_depth and int(s_depth) < int(min_depth):
                        continue
                    if max_depth and int(s_depth) > int(max_depth):
                        continue
                ssr_type = item[7].split("ST=")[1]  # 计算SSR类型分布
                if item[0] not in stat_info[cs].keys():
                    stat_info[cs][item[0]] = {"ssr": 0, "c": 0, "c*": 0, "p1": 0, "p2": 0, "p3": 0, "p4": 0, "p5": 0, "p6": 0}
                stat_info[cs][item[0]]["ssr"] += 1
                stat_info[cs][item[0]][ssr_type] += 1
                file_info[cs]["detail"].write("\t".join(item[:6]) + "\t" + str(c_count) + "\t" + item[sample_index[c_sample]].split(":")[-1])
                file_info[cs]["detail"].write("\t" + str(s_count) + "\t" + item[sample_index[s_sample]].split(":")[-1] + "\n")
    for cs in sample_info:
        for chr in stat_info[cs].keys():
            file_info[cs]["stat"].write(chr + "\t" + str(stat_info[cs][chr]["ssr"]) + "\t" + str(stat_info[cs][chr]["c"]) + "\t"
                                          + str(stat_info[cs][chr]["c*"]) + "\t" + str(stat_info[cs][chr]["p1"]) + "\t" + str(stat_info[cs][chr]["p2"])
                                          + "\t" + str(stat_info[cs][chr]["p3"]) + "\t" + str(stat_info[cs][chr]["p4"]) + "\t"
                                          + str(stat_info[cs][chr]["p5"]) + "\t" + str(stat_info[cs][chr]["p6"]) + "\n")


def multiple_ssr_compare(ssr_file, config_file, output_dir):
    """
    多条件样本间SSR比较
    params: {"region_select": {"chr1": 1-10000}, "allele_num": 2/3, "sample_params": {"sample1_VS_sample2":
             {"sample1":{"gene_type":"same/mixed","depth":"-100"}, "sample2":{"gene_type":"same/mixed","depth":"-100"},
             "compare_type": "same/diff"}}, "group_params": {"group1": {"sample":["sample1", "sample2"],"depth":,"gene_fre":-,"miss_fre":-}}}
            allele_num: 等位基因数目(2/3)
            compare_type: 比较类型(same/diff)
            gene_type: 基因型(same/mixed)
            depth: 测序深度(-分隔)
            sample_params: 样本比较的参数
            group_params: 样本组比较的参数
    """
    params = json.loads(open(config_file).read())
    region_select = params["region_select"]
    sample_info = params["sample_params"] if "sample_params" in params.keys() else None
    group_info = params["group_params"] if "group_params" in params.keys() else None
    stat_info, sample_index = {}, {}
    ssr_detail = open(os.path.join(output_dir, "ssr_detail.xls"), "w")
    ssr_stat = open(os.path.join(output_dir, "ssr_stat.xls"), "w")
    ssr_detail.write("#CHROM\tStart\tEnd\tSSR ID\tRepeat unit\tReference Repeat count")
    if sample_info:
        sample_list = []
        for cs in sample_info.keys():
            sample_list.append(cs.split("_VS_")[0])
            sample_list.append(cs.split("_VS_")[1])
        sample_list = list(set(sample_list))
        for sample in sample_list:
            ssr_detail.write("\t" + sample + "_repeat_count\t" + sample + "_allele_depth")
    if group_info:
        group_list = []
        for group in group_info.keys():
            group_list.append(group)
            ssr_detail.write("\t" + group + "_allele_frequency")
    ssr_detail.write("\n")
    ssr_stat.write("#CHROM\tSSR Number\tc\tc*\tp1\tp2\tp3\tp4\tp5\tp6\n")
    with open(ssr_file, "rb") as f:
        for line in f:
            if line.startswith("##"):
                continue
            item = line.strip().split("\t")
            if line.startswith("#CHROM"):  # 样本的index
                for i in range(9, len(item)):
                    sample_index[item[i]] = i
            if item[6] != "PASS":
                continue
            if region_select != "all":  # 区域过滤
                status = region_filter(item, region_select)
                if not status:
                    continue
            # ref_gene_type = list(set(re.findall("\d+", item[4])))
            # if len(ref_gene_type) not in allele_num:  # 等位基因数目过滤(过滤参考基因组)
            #     continue
            if sample_info:
                status, sample_stat = sample_filter(sample_info, sample_index, item)
                if not status:
                    continue
            if group_info:
                status, group_stat = group_filter(group_info, sample_index, item)
                if not status:
                    continue
            ssr_type = item[7].split("ST=")[1]  # 计算SSR类型分布
            if item[0] not in stat_info.keys():
                stat_info[item[0]] = {"ssr": 0, "c": 0, "c*": 0, "p1": 0, "p2": 0, "p3": 0, "p4": 0, "p5": 0, "p6": 0}
            stat_info[item[0]]["ssr"] += 1
            stat_info[item[0]][ssr_type] += 1
            ssr_detail.write("\t".join(item[:6]))
            if sample_info:
                for s in sample_list:
                    ssr_detail.write("\t" + str(sample_stat[s][0]) + "\t" + str(sample_stat[s][1]))
            if group_info:
                for g in group_list:
                    ssr_detail.write("\t" + str(group_stat[g][1]))
            ssr_detail.write("\n")
    for chr in stat_info.keys():
        ssr_stat.write(chr + "\t" + str(stat_info[chr]["ssr"]) + "\t" + str(stat_info[chr]["c"]) + "\t"
                       + str(stat_info[chr]["c*"]) + "\t" + str(stat_info[chr]["p1"]) + "\t" + str(stat_info[chr]["p2"])
                       + "\t" + str(stat_info[chr]["p3"]) + "\t" + str(stat_info[chr]["p4"]) + "\t"
                       + str(stat_info[chr]["p5"]) + "\t" + str(stat_info[chr]["p6"]) + "\n")


parser = argparse.ArgumentParser(description="SSR比较分析")
parser.add_argument("-ssr", "--ssr_file", required=True)
parser.add_argument("-c", "--config_file", required=True)
parser.add_argument('-m', '--model', help='single or multiple', required=True)
parser.add_argument("-o", "--output_dir", required=True)

args = vars(parser.parse_args())

if args["model"] == "single":
    single_ssr_compare(args["ssr_file"], args["config_file"], args["output_dir"])
elif args["model"] == "multiple":
    multiple_ssr_compare(args["ssr_file"], args["config_file"], args["output_dir"])
else:
    raise Exception("参数-m类型只能是single或multipe")


"""
单条件批量参数: {"allele_num": 2, "compare_type": "same", "sample_info": ["BDZ_VS_AA"], "region_select": {"chr1": "1-500000"}, "depth": "2-100", "gene_type": "same"}
多条件组合参数：{"allele_num": 2, "sample_params": {"BDZ_VS_AA": {"AA": {"depth": "-100", "gene_type": "same"}, "compare_type": "diff", "BDZ": {"depth": "-100", "gene_type": "same"}}}, "group_params": {"group1": {"sample": ["BDZ", "AA"], "depth": "-", "gene_fre": "-0.7", "miss_fre": "-"}}, "region_select": "all"}
"""
