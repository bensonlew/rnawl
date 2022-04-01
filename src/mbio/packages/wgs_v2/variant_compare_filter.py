# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 2019.02.21
import re
import argparse
import time
from collections import defaultdict

parse = argparse.ArgumentParser(description="用于变异位点比较分析接口的package")
parse.add_argument("-i", "--vcf", help="输入用于分析的vcf文件", required=True)
parse.add_argument("-c", "--config", help="输入config文件", required=True)
parse.add_argument("-o", "--output", help="输出结果目录", required=True)
parse.add_argument("-g", "--group", help="输出结果目录", required=True)
parse.add_argument("-n", "--name", help="生成文件名称前缀", required=True)
args = vars(parse.parse_args())
vcf_path = args["vcf"]
pop_filter_path = args["output"] + "/" + args["name"] + ".filter.vcf"
config_path = args["config"]
group_list = args["group"]
eff_path = args["output"] + "/" + args["name"] + ".eff"
func_path = args["output"] + "/" + args["name"] + ".func"
detail_path = args["output"] + "/" + args["name"] + ".detail"


def sort_key(list1):
    list_chr_d = []
    list_sca_d = []
    list_chr_D = []
    list_sca_D = []
    for chr in list1:
        if chr.startswith("chr"):
            if len(re.findall("\d+", chr)) != 0:
                list_chr_d.append(chr)
            else:
                list_chr_D.append(chr)
        elif chr.startswith("sca"):
            if len(re.findall("\d+", chr)) != 0:
                list_sca_d.append(chr)
            else:
                list_sca_D.append(chr)
        else:
            raise Exception("染色体中不合规范，请核实！")
    list_chr_d.sort(key=lambda k: int(re.findall('\d+', k)[0]))
    list_sca_d.sort(key=lambda k: int(re.findall('\d+', k)[0]))
    list_chr_D.sort()
    list_sca_D.sort()
    new_list = list_chr_d + list_chr_D + list_sca_d + list_sca_D
    return new_list


def check_chr(chr1, pos1, chr1_list, chr_list1_new, region_list1):
    if not chr1_list:
        return 2
    else:
        if chr1 != chr_list1_new[-1] and sort_key([chr1, chr_list1_new[-1]])[1] == chr1:
            return 0
        else:
            for chr in chr1_list:
                if chr1 == chr:
                    range_new_min = region_list1[chr1_list.index(chr)].strip().split(":")[1].split("-")[0]
                    range_new_max = region_list1[chr1_list.index(chr)].strip().split(":")[1].split("-")[1]
                    if int(pos1) > int(range_new_max):
                        return 1
                    elif int(range_new_max) >= int(pos1) >= int(range_new_min):
                        return 2
                    else:
                        return 1
                else:
                    if chr == chr_list1_new[-1]:
                        return 1


def variant_filter(variant_type1, allele_num1, variant_eff, gt1, gt2, dep1, dep2, dep3, dep4, genotype1, ad1, ad2,
                   max_miss, maf1, maf2, variant_type, allele_num, eff, locations, gt_dict, genotype_dict, dp_dict,
                   alle_dp_dict, group_dict, row, sample_list1, group_list1, region_list1, gts, deps,
                   sample_compare_list1):
    """
    第一过滤的函数，vcf中的每一行都需要来这里过滤，通过的执行下一步，不能通过的停止。
    :return:
    """
    n = 1  # 用于最终判断是否通过测试
    check = 'no'  # 用于判断能否通过region的测试
    group_maf = {}  # 记录组的平均基因型频率
    if variant_type not in variant_type1:  # 判断snp,indel 是否符合条件
        return False
    else:
        if str(allele_num) not in allele_num1:  # 判断等位基因数目是否符合条件
            return False
        else:
            if len(list(set(variant_eff) & set(eff))) == 0:  # 求两个集合的交集，如果有交集说明可以通过考核。
                return False
            else:
                # if genome_region == "all":    # 这里还需要补充不同的region的格式。
                #     pass
                for region in region_list1:
                    if region == "all":  # 当传入为all的时候相当于本条命令不存在。
                        check = "yes"
                        break
                    else:
                        if not locations.strip().split("_")[0] in region.strip().split(":")[0]:  # 染色体不在范围内
                            pass
                        else:
                            if not int(locations.strip().split("_")[1]) in range(
                                    int(region.strip().split(":")[1].split("-")[0]), int(
                                        region.strip().split(":")[1].split("-")[1]) + 1):  # 基因位置不在范围内, 必须得+1
                                pass
                            else:
                                check = "yes"
                if check == "yes":
                    # print dp_dict
                    if len(sample_list1) != 0:
                        for samples in sample_compare_list1:
                            sample1, sample2 = samples.split('|')
                            if sample1.lower() != 'reference':
                                if sample_name1.lower() == "reference":
                                    pass
                                else:
                                    if gt_dict[sample_name1] == "null":
                                        return False
                            if sample2.lower() != 'reference':
                                if sample_name2.lower() == "reference":
                                    pass
                                else:
                                    if gt_dict[sample_name2] == "null":
                                        return False
                            if sample1.lower() == 'reference':
                                if sample2.lower() == 'reference':
                                    raise Exception('两个样本不能同时为reference')
                                if gt_dict[sample2] == gts[sample2] or gts[sample2] == 'homo, hete':  # 判断纯合与杂合
                                    dep1, dep2 = deps[sample2].split(',')
                                    try:
                                        if int(dep1) <= int(dp_dict[sample2]) <= int(dep2):
                                            if genotype1[samples] == 'all':
                                                n = 0
                                            elif genotype1[samples] == 'same':
                                                if int(genotype_dict[sample2][0]) == 0:
                                                    n = 0
                                            elif genotype1[samples] == 'diff':
                                                if int(genotype_dict[sample2][0]) == 1:
                                                    n = 0
                                    except:
                                        return False
                            elif sample2.lower() == 'reference':
                                if sample1.lower() == 'reference':
                                    raise Exception('两个样本不能同时为reference')
                                if gt_dict[sample1] == gts[sample1] or gts[sample1] == 'homo, hete':  # 判断纯合与杂合
                                    dep1, dep2 = deps[sample1].split(',')
                                    try:
                                        if int(dep1) <= int(dp_dict[sample1]) <= int(dep2):
                                            if genotype1[samples] == 'all':
                                                n = 0
                                            elif genotype1[samples] == 'same':
                                                if int(genotype_dict[sample1][0]) == 0:
                                                    n = 0
                                            elif genotype1[samples] == 'diff':
                                                if int(genotype_dict[sample1][0]) == 1:
                                                    n = 0
                                    except:
                                        return False
                            else:
                                if gt_dict[sample1] == 'null' or gt_dict[sample2] == 'null':
                                    return False
                                else:
                                    if gt_dict[sample1] == gts[sample1] and gt_dict[sample2] == gts[sample2]:
                                        # 纯合杂合all正确
                                        dep1, dep2 = deps[sample1].split(',')
                                        dep3, dep4 = deps[sample2].split(',')
                                        try:
                                            if int(dep1) <= int(dp_dict[sample1]) <= int(dep2) and int(dep3) <= \
                                                    int(dp_dict[sample2]) <= int(dep4):  # 这里如果是缺失没有值的话怎么么处理，上面已经做了处理。
                                                if genotype1[samples] == "all":
                                                    n = 0  # 当n赋值为0的时候，意味着返回值为True
                                                elif genotype1[samples] == "same":
                                                    if set(genotype_dict[sample1]) == set(genotype_dict[sample2]):
                                                        n = 0  # 当n赋值为0的时候，意味着返回值为True
                                                elif genotype1[samples] == "diff":
                                                    if set(genotype_dict[sample1]) != set(genotype_dict[sample2]):
                                                        n = 0  # 当n赋值为0的时候，意味着返回值为True
                                        except:
                                            return False
                                    elif (gts[sample1] == 'homo, hete' and gt_dict[sample2] == gts[sample2]) \
                                            or (
                                            gts[sample2] == 'homo, hete' and gt_dict[sample1] == gts[sample1]) or \
                                            (gts[sample1] == 'homo, hete' and gts[sample2] == 'homo, hete'):
                                        # add by hd 解决基因型为all的时候bug  临时添加后面彬彬再优化
                                        dep1, dep2 = deps[sample1].split(',')
                                        dep3, dep4 = deps[sample2].split(',')
                                        try:
                                            if int(dep1) <= int(dp_dict[sample1]) <= int(dep2) and int(dep3) <= \
                                                    int(dp_dict[sample2]) <= int(dep4):  # 这里如果是缺失没有值的话怎么么处理，上面已经做了处理。
                                                if genotype1[samples] == "all":
                                                    n = 0  # 当n赋值为0的时候，意味着返回值为True
                                                elif genotype1[samples] == "same":
                                                    if set(genotype_dict[sample1]) == set(genotype_dict[sample2]):
                                                        n = 0  # 当n赋值为0的时候，意味着返回值为True
                                                elif genotype1[samples] == "diff":
                                                    if set(genotype_dict[sample1]) != set(genotype_dict[sample2]):
                                                        n = 0  # 当n赋值为0的时候，意味着返回值为True
                                        except:
                                            return False

                    # if len(sample_list1) != 0:
                    #     if gt_dict[sample_name1] == "null" or gt_dict[sample_name2] == "null":  # 排除没有检测出的位点。
                    #         return False
                    #     else:
                    #         if gt_dict[sample_name1] in gt1 and gt_dict[sample_name2] in gt2:  # 纯合杂合all正确
                    #             if int(dep1) <= int(dp_dict[sample_name1]) <= int(dep2) and int(dep3) <= \
                    #                     int(dp_dict[sample_name2]) <= int(dep4):  # 这里如果是缺失没有值的话怎么么处理，上面已经做了处理。
                    #                 if genotype1 == "all":
                    #                     n = 0  # 当n赋值为0的时候，意味着返回值为True
                    #                 elif genotype1 == "same":
                    #                     if set(genotype_dict[sample_name1]) == set(genotype_dict[sample_name2]):
                    #                         n = 0  # 当n赋值为0的时候，意味着返回值为True
                    #                 elif genotype1 == "diff":
                    #                     if set(genotype_dict[sample_name1]) != set(genotype_dict[sample_name2]):
                    #                         n = 0  # 当n赋值为0的时候，意味着返回值为True
                    if len(group_list1) != 0:
                        n = 1
                        for keys in group_dict.keys():
                            miss_value = 0
                            ex = 0
                            sample_list_new = []  # 用于存储通过过滤的样本名，然后这里的名称用于计算基因频率。
                            for samples in group_dict[keys]:
                                if samples not in dp_dict.keys():  # 如果不存在说明此处存在缺失的情况。
                                    miss_value += 1
                                else:
                                    if not int(ad1) <= int(dp_dict[samples]) <= int(ad2):  # 深度值不符合，这里要做缺失率处理。
                                        miss_value += 1
                                    else:
                                        ex += 1
                                        sample_list_new.append(samples)
                            miss_ratio = float(miss_value) / float(len(group_dict[keys]))
                            if float(miss_ratio) <= float(max_miss):  # 按照现在的逻辑，缺失率就同时吧深度也过滤掉了
                                gene_num = 0  # 用于储存基因型是0的个数。
                                if not len(sample_list_new) == 0:
                                    for i in sample_list_new:
                                        for gene in genotype_dict[i]:
                                            if gene == "0":
                                                gene_num += 1
                                    global gene_fre
                                    gene_fre = float(gene_num) / (float(len(sample_list_new)) * 2)
                                    if float(maf1) <= gene_fre <= float(maf2):
                                        group_maf[keys] = gene_fre
                                        n = 0
                                    else:
                                        n = 1
                                        break
                                else:
                                    n = 1
                                    break
                            else:
                                n = 1
                                break
                    if n == 0:
                        write_detail(row, genotype_dict, alle_dp_dict, sample_list1, group_list1, group_maf,
                                     variant_type)
                        return True
                    else:
                        return False


def write_vcf(row):
    with open(pop_filter_path, "a+")as b:
        b.write(row)


def write_detail(row, genotype_dict, alle_dp_dict, sample_list1, group_list1, group_maf,
                 variant_type):  # 第二个，第三个可以随便传给我，我需要知道名字
    with open(detail_path, "a+") as w2:
        split1 = row.strip().split("\t")
        row_new = split1[0:5]
        row_new.append(variant_type)
        anno = [split1[7].strip().split(";")[3].strip().split(",")[m].strip().split("=")[1] if m == 0 else
                split1[7].strip().split(";")[3].strip().split(",")[m]
                for m in range(len(split1[7].strip().split(";")[3].strip().split(",")))]
        row_new.append(";".join(anno))
        """
        pos用于存储位置信息供前端取用。因为indel可能会有多个，这里以第一个进行计算"""
        pos = str(split1[0]) + ":" + str(split1[1]) + "-" + str(
            int(split1[1]) + len(split1[4].strip().split(",")[0]) -
            1)
        row_new.append(pos)
        if len(sample_list1) == 0:
            pass
        else:
            for i in sample_list1:
                if i.lower() == 'reference':
                    pass
                else:
                    new_list = []  # 用于存放新的基因型
                    for j in genotype_dict[i]:
                        if j == "0":
                            new_list.append(split1[3])
                        else:
                            new_list.append(split1[4].split(",")[int(j) - 1])
                    row_new.append("/".join(new_list))
                    row_new.append(alle_dp_dict[i])
        if len(group_list1) == 0:
            pass
        else:
            for j in group_list1:
                row_new.append(str(group_maf[j]))
        w2.write("\t".join(row_new) + "\n")


def write_fun_eff(func):  # 功能和功效统一在这里计算。

    for func_item in set(func):  # 把list转为话set的话，可以防止重复计数
        if variant_type == "SNP":
            if func_item not in dict_snp_func.keys():
                dict_snp_func[func_item] = 1
            else:
                dict_snp_func[func_item] += 1
        else:
            if func_item not in dict_indel_func.keys():
                dict_indel_func[func_item] = 1
            else:
                dict_indel_func[func_item] += 1
    with open(func_path, "w") as w1:
        w1.write("effect type\tsnp number\tindel number\ttotal number\n")
        list_key = []  # 建立一个新的list用于存储indel和snp两个的key数值。
        for key1 in dict_snp_func.keys():
            if key1 not in list_key:
                list_key.append(key1)
        for key2 in dict_snp_func.keys():
            if key2 not in list_key:
                list_key.append(key2)
        for item in list_key:
            if item not in dict_snp_func.keys():
                dict_snp_func[item] = 0
            if item not in dict_indel_func.keys():
                dict_indel_func[item] = 0
            w1.write(item + "\t" + str(dict_snp_func[item]) + "\t" + str(dict_indel_func[item]) + "\t" +
                     str(dict_snp_func[item] + dict_indel_func[item]) + "\n")


def write_eff(snp_high, indel_high, snp_moderate, indel_moderate, snp_low, indel_low, snp_modifier, indel_modifier):
    with open(eff_path, "w") as w:
        w.write("impact type\tsnp number\tindel number\ttotal number\n")
        w.write("high\t" + str(snp_high) + "\t" + str(indel_high) + "\t" + str(snp_high + indel_high) + "\n")
        w.write(
            "moderate\t" + str(snp_moderate) + "\t" + str(indel_moderate) + "\t" + str(snp_moderate + indel_moderate
                                                                                       ) + "\n")
        w.write("low\t" + str(snp_low) + "\t" + str(indel_low) + "\t" + str(snp_low + indel_low) + "\n")
        w.write(
            "modifier\t" + str(snp_modifier) + "\t" + str(indel_modifier) + "\t" + str(snp_modifier + indel_modifier
                                                                                       ) + "\n")


snp_high = 0
snp_moderate = 0
snp_low = 0
snp_modifier = 0
indel_high = 0
indel_moderate = 0
indel_low = 0
indel_modifier = 0
dict_snp_func = {}  # 用于储存func的元素
dict_indel_func = {}
dep1 = 0
dep2 = 0
dep3 = 0
dep4 = 0
gt1 = ""
gt2 = ""
genotype1 = {}
ad1 = 0
ad2 = 0
max_miss = 0
maf1 = 0
maf2 = 0
sample_list1 = []
group_list1 = []
region_list1 = []
chr_list1 = []
deps = {}
gts = {}
sample_compare_list1 = []

"""首先解析配置文件
"""
with open(config_path)as f:
    n = 1
    lines = f.readlines()
    variant_type1 = lines[0].strip().split("=")[1].strip().split(",")
    allele_num1 = lines[2].strip().split("=")[1].strip().split(",")
    variant_eff = lines[1].strip().split("=")[1].strip().split(",")
    # genome_region = lines[3].strip("=")[1]
    for line in lines[3:]:
        if re.match("Region", line):
            genome_region = line.strip().split("=")[1]
            if genome_region not in region_list1:
                region_list1.append(genome_region)
                if genome_region != "all":
                    chr_region = genome_region.strip().split(":")[0]
                    chr_list1.append(chr_region)
        if re.match("Sample Diff", line):
            sample_name1 = line.strip().split("=")[1].split(",")[0]
            sample_name2 = line.strip().split("=")[1].split(",")[4]
            if sample_name1 not in sample_list1:
                sample_list1.append(sample_name1)
            if sample_name2 not in sample_list1:
                sample_list1.append(sample_name2)
            sample_compare_list1.append(sample_name1 + '|' + sample_name2)
            dep1 = line.strip().split("=")[1].split(",")[1]
            dep2 = line.strip().split("=")[1].split(",")[2]
            deps[sample_name1] = ','.join([dep1, dep2])
            dep3 = line.strip().split("=")[1].split(",")[5]
            dep4 = line.strip().split("=")[1].split(",")[6]
            deps[sample_name2] = ','.join([dep3, dep4])
            gt1 = "homo, hete" if line.strip().split("=")[1].split(",")[3] == "all" else \
                line.strip().split("=")[1].split(",")[3]  # 表示是homo, heter, 或者是all
            gt2 = "homo, hete" if line.strip().split("=")[1].split(",")[7] == "all" else \
                line.strip().split("=")[1].split(",")[7]
            gts[sample_name1] = gt1
            gts[sample_name2] = gt2
            genotype1[sample_name1 + '|' + sample_name2] = line.strip().split("=")[1].split(",")[8]  # 表示是相同还是不相同
        if re.match("Group Info", line):
            group_name = line.strip().split("=")[1].split(",")[0]
            group_list1.append(group_name)
            ad1 = line.strip().split("=")[1].split(",")[1]
            ad2 = line.strip().split("=")[1].split(",")[2]
            max_miss = line.strip().split("=")[1].split(",")[3]
            maf1 = line.strip().split("=")[1].split(",")[4]
            maf2 = line.strip().split("=")[1].split(",")[5]

sort_key(chr_list1)
header_list = ["#chr", "pos", "snp/indel_id", "tref_geno", "alt", "type", "annotation", "pos"]
# 后面一个pos是将染色体和开始结束位置连在一起的一个pos
hasReference = False
for i in sample_list1:
    if i.lower() == 'reference':
        hasReference = True
        continue
    header_list.append(i+"_genotype")
    header_list.append(i + "_allele_depth")
for j in group_list1:
    header_list.append(j+"_alle_fre")
with open(detail_path, "a+") as w2:
    for i in header_list:
        w2.write(i+"\t")
    w2.write("\n")
"""
打开vcf文件进行一行行分析"""
with open(vcf_path) as m:
    rows = m.readlines()
    sample_list = []
    i = 0
    for row in rows:
        i += 1
        gt_dict = {}  # 用于存储每一行中每一个样本的基因型，homo，heter等
        dp_dict = {}  # 用于存放每一行中每一个样本的深度信息
        genotype_dict = {}  # 用于存储每一行中的每一个样本的基因型信息，atcg等
        genotype = []
        alle_dp_dict = {}  # 用于存放等位基因深度信息，详情表中需要使用。
        if re.match("##", row):
            write_vcf(row)
        elif re.match("#CHROM", row):
            write_vcf(row)
            sample_list = row.strip().split("\t")[9:]
        else:
            variant = row.strip().split("\t")
            chr1 = variant[0]
            pos1 = variant[1]
            index = 0  # 直接断掉
            return_code = check_chr(chr1, pos1, chr_list1, sort_key(chr_list1), region_list1)
            if return_code == 0:
                break
            elif return_code == 1:
                continue
            else:
                gene_list = (str(variant[3]) + "," + str(variant[4])).strip().split(",")
                variant_type = "SNP"
                for gene in gene_list:
                    if len(gene) > 1:
                        variant_type = "INDEL"
                        break
                    else:
                        pass
                allele_num = len(variant[4].strip().split(","))+1
                # modified by Binbin Zhao@20190910, eff,func在第七列中的位置并不固定。
                n = 0
                eff_func = variant[7].strip().split(";")
                for k in range(len(eff_func)):
                    if re.match("ANN=", eff_func[k]):
                        n = k
                        break
                eff = [x.strip().split("|")[2] for x in variant[7].strip().split(";")[n].strip().split(",")]
                func = [x.strip().split("|")[1] for x in variant[7].strip().split(";")[n].strip().split(",")]
                locations = variant[0]+"_" + variant[1]
                format1 = variant[8].strip().split(":")
                genotype_origin = [y.strip().split(":")[format1.index("GT")] for y in variant[9:]]
                genotype_final = [x.strip().split("/") for x in genotype_origin]
                try:
                    depth = [y.strip().split(":")[format1.index("DP")] for y in variant[9:]]
                    alle_depth = [y.strip().split(":")[format1.index("AD")] for y in variant[9:]]
                except:
                    pass
                for k in genotype_origin:
                    if k == "./." or k == ".":
                        k1 = "null"
                        genotype.append(k1)
                    elif len(set(k.strip().split("/"))) == 1:
                        k1 = "homo"
                        genotype.append(k1)
                    elif len(set(k.strip().split("/"))) > 1:
                        k1 = "hete"
                        genotype.append(k1)
                    else:
                        raise Exception("genotype未知错误，请核实！")
                for (x, y) in zip(sample_list, genotype):
                    gt_dict[x] = y                                 # 用来储存纯合还是杂合
                for (y, z) in zip(sample_list, depth):
                    if not gt_dict[y] == "null":
                        dp_dict[y] = z                            # 用于储存深度信息
                for (x, y) in zip(sample_list, genotype_final):
                    if not gt_dict[x] == "null":
                        genotype_dict[x] = y                      # 用于存储基因型，内容为0,1
                for (m, n) in zip(sample_list, alle_depth):
                    if not gt_dict[m] == "null":
                        alle_dp_dict[m] = n  # 用于存储基因型，内容为0,1
                # print dp_dict
                with open(group_list) as g:
                    group_dict = {}
                    group_info = g.readlines()
                    for k2 in group_info:
                        group_name1 = k2.strip().split(":")[0]
                        group_sample = k2.strip().split(":")[1].strip().split(",")
                        if group_name1 in group_list1:
                            if group_name1 not in group_dict.keys():
                                group_dict[group_name1] = group_sample
                if variant_filter(variant_type1, allele_num1, variant_eff, gt1, gt2, dep1, dep2, dep3, dep4, genotype1,
                                  ad1, ad2, max_miss, maf1, maf2, variant_type, allele_num, eff, locations, gt_dict,
                                  genotype_dict, dp_dict, alle_dp_dict, group_dict, row, sample_list1, group_list1,
                                  region_list1, gts, deps, sample_compare_list1):
                    if "HIGH" in eff:
                        if variant_type == "SNP":
                            snp_high += 1
                        else:
                            indel_high += 1
                    if "MODERATE" in eff:
                        if variant_type == "SNP":
                            snp_moderate += 1
                        else:
                            indel_moderate += 1
                    if "LOW" in eff:
                        if variant_type == "SNP":
                            snp_low += 1
                        else:
                            indel_low += 1
                    if "MODIFIER" in eff:
                        if variant_type == "SNP":
                            snp_modifier += 1
                        else:
                            indel_modifier += 1
                    write_vcf(row)
                    write_eff(snp_high, indel_high, snp_moderate, indel_moderate, snp_low, indel_low, snp_modifier,
                              indel_modifier)
                    write_fun_eff(func)
                else:
                    pass

data_dict = {}  # 保存不同chr对应的格式值{chr1: [1,2]}, 第一个值代表snp个数，第二值代表indel个数
snp_indeloutfile = args["output"] + "/" + args["name"] + ".snp_indel_stat.txt"
with open(pop_filter_path, 'r') as r, open(snp_indeloutfile, 'w') as w:
    w.write('#chr id\tsnp num\t indel num\t total num\n')
    for line in r:
        if re.match('#', line):
            pass
        else:
            temp = line.strip().split('\t')
            gene_list = (str(temp[3]) + "," + str(temp[4])).strip().split(",")
            variant_type = "SNP"
            for gene in gene_list:
                if len(gene) > 1:
                    variant_type = "INDEL"
                    break
            if temp[0] not in data_dict.keys():
                if variant_type == "SNP":
                    data_dict[temp[0]] = [1, 0]
                else:
                    data_dict[temp[0]] = [0, 1]
            else:
                if variant_type == "SNP":
                    data_dict[temp[0]][0] += 1
                else:
                    data_dict[temp[0]][1] += 1
    for keys in sorted(data_dict.keys()):
        w.write('{}\t{}\t{}\t{}\n'.format(keys, data_dict[keys][0], data_dict[keys][1], sum(data_dict[keys])))


