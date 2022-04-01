# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# lasted modified by hongdong@20190314

import re
import argparse


def consensus_stat(input, stat_path, cover_path, pop_tag, total_sample):
    """
    统计consensus的结果
    lasted modified by hongdong@20190314
    1)所有的tab转为空格
    2）Consensus聚类覆盖度分布根据样本个数进行动态结果展示
    """
    total_consensus, total_leng, total_depth = 0, 0, 0
    cover_20, cover_50, cover_80 = 0, 0, 0
    sample_20 = total_sample * 0.2
    sample_50 = total_sample * 0.5
    sample_80 = total_sample * 0.8
    sample_num = 0
    sample_num_list = {}
    seq_list = []
    with open(input, "r") as f, open(stat_path, "w") as w, open(cover_path, "w") as w1, open(pop_tag, "w") as w2:
        for line in f:
            if line.startswith("//"):
                if sample_num not in sample_num_list.keys():
                    sample_num_list[sample_num] = 0
                sample_num_list[sample_num] += 1
                if sample_num >= sample_20:
                    cover_20 += 1
                if sample_num >= sample_50:
                    cover_50 += 1
                if sample_num >= sample_80:
                    cover_80 += 1
                m = re.match(r"//.+\|(\d+)\|", line)
                num = len(seq_list[0])
                seq, snp = "", ""
                for i in range(num - leng -1, num):
                    seq += seq_list[0][i]
                    if line[i] in ["*", "-"]:
                        snp += "*"
                    else:
                        snp += " "
                w2.write(">locus_" + m.group(1) + "\n" + seq + str(snp) + "\n")
                total_consensus += 1
                total_leng += leng
                sample_num = 0
                seq_list = []
            else:
                total_depth += 1
                sample_num += 1
                l = re.match(r"\S+\s+(\S+)", line)
                leng = len(l.group(1))
                seq_list.append(line)
        # average_leng = round(float(total_leng) / total_consensus, 4)
        average_depth = round(float(total_depth) / total_consensus / total_sample, 4) * 100
        w.write("Consensus Number\t{}\n".format(str(total_consensus)))
        # w.write("Average Length\t{}\n".format(str(average_leng)))
        # w.write("Average Depth\t{}\n".format(str(average_depth)))
        w.write("Average Coverage\t{}%\n".format(str(average_depth)))
        w.write("Consensus Coverage (20%)\t{}\n".format(str(cover_20)))
        w.write("Consensus Coverage (50%)\t{}\n".format(str(cover_50)))
        w.write("Consensus Coverage (80%)\t{}\n".format(str(cover_80)))
        sample_num_list_ = sample_num_list.keys()
        w1.write("#Coverage Sample Num\tConsensus Number\n")
        if int(total_sample) <= 10:
            for num in sorted(sample_num_list_):
                rate = float(num) / float(total_sample) * 100
                w1.write(str(round(rate, 2)) + "%\t" + str(sample_num_list[num]) + "\n")
        else:
            cov_0_10, cov_10_20, cov_20_30, cov_30_40, cov_40_50, cov_50_60, cov_60_70, cov_70_80, cov_80_90, cov_90_100 \
                = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            for num in sorted(sample_num_list_):
                rate = float(num) / float(total_sample) * 100
                if 0 < rate < 10:
                    cov_0_10 += int(sample_num_list[num])
                elif 10 <= rate < 20:
                    cov_10_20 += int(sample_num_list[num])
                elif 20 <= rate < 30:
                    cov_20_30 += int(sample_num_list[num])
                elif 30 <= rate < 40:
                    cov_30_40 += int(sample_num_list[num])
                elif 40 <= rate < 50:
                    cov_40_50 += int(sample_num_list[num])
                elif 50 <= rate < 60:
                    cov_50_60 += int(sample_num_list[num])
                elif 60 <= rate < 70:
                    cov_60_70 += int(sample_num_list[num])
                elif 70 <= rate < 80:
                    cov_70_80 += int(sample_num_list[num])
                elif 80 <= rate < 90:
                    cov_80_90 += int(sample_num_list[num])
                elif 90 <= rate <= 100:
                    cov_90_100 += int(sample_num_list[num])
            w1.write("0%-10%\t{}\n".format(cov_0_10))
            w1.write("10%-20%\t{}\n".format(cov_10_20))
            w1.write("20%-30%\t{}\n".format(cov_20_30))
            w1.write("30%-40%\t{}\n".format(cov_30_40))
            w1.write("40%-50%\t{}\n".format(cov_40_50))
            w1.write("50%-60%\t{}\n".format(cov_50_60))
            w1.write("60%-70%\t{}\n".format(cov_60_70))
            w1.write("70%-80%\t{}\n".format(cov_70_80))
            w1.write("80%-90%\t{}\n".format(cov_80_90))
            w1.write("90%-100%\t{}\n".format(cov_90_100))


parser = argparse.ArgumentParser(description='统计consensus的信息')
parser.add_argument('-i', '--input', help='data.loci', required=True)
parser.add_argument('-s', '--stat_path', help='output consensus_stat.xls', required=True)
parser.add_argument('-c', '--cover_path', help='output consensus_coverage.xls', required=True)
parser.add_argument('-p', '--pop_tag', help='output populations.tag', required=True)
parser.add_argument('-t', '--total_sample', help='total sample num', required=True)

args = vars(parser.parse_args())
consensus_stat(args['input'], args['stat_path'], args['cover_path'], args['pop_tag'], int(args['total_sample']))
