## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "HONGDONG"
# last_modify:20190313

import re
import argparse
import os
import time
from collections import defaultdict


class WindowSliding(object):
    def __init__(self, infile, step, output_dir):
        """
        :param infile: 位点信息文件
        :param step: 划窗的步长
        :param output_dir: 运行结果目录
        """
        self.infile = infile
        self.step = step
        self.output_dir = output_dir

    def sliding(self):
        """
        算法逻辑：
        1）获取要进行划窗的列表的最大值，使用range，按照step将列表分成需求等分。
        2）然后将列表的每个字节与等分进行比较大小，这个时候生成一个列表【10， 20， 30】，这里20，其实是包含了10的结果了。
        3）shuzuzhuanhua将数据中后面一个数减去前面一个数，就得到我们想要的结果了。
        :return:
        """
        win_data = {}
        chr_data = self.read_vcf_file()
        outfile = os.path.join(self.output_dir, "win_{}_result.txt".format(self.step))
        with open(outfile, 'w') as w:
            w.write('#chrname\twin_result\n')
            for key in sorted(chr_data.keys()):
                pos = chr_data[key]
                win_data[key] = defaultdict(list)
                pos.sort()
                end = pos[-1]
                win_result = []
                regions = range(0, end + self.step, self.step)[1:]
                for n in pos:
                    for i in range(0, len(regions)):
                        if regions[i] not in win_data[key].keys():
                            win_data[key][regions[i]] = 0
                        if n <= regions[i]:
                            win_data[key][regions[i]] += 1
                            continue
                for keys in sorted(win_data[key].keys()):
                    if keys != 0:
                        win_result.append(win_data[key][keys])
                # print self.shuzuzhuanhua(win_result)
                w.write("{}\t{}\n".format(key, ','.join(self.shuzuzhuanhua(win_result))))

    def shuzuzhuanhua(self, list_):
        """
        将数据中后面一个数减去前面一个数
        :param list_:
        :return:
        """
        new_list = []
        for i in range(0, len(list_)):
            if i == 0:
                new_list.append(str(list_[0]))
            else:
                new_list.append(str(list_[i] - list_[i - 1]))
        return new_list

    def read_vcf_file(self):
        """
        文件按照{'chr1': [1000,2000,3000], 'chr2': [2000,3000,5000]}存放
        输入的文件可以是chr start（snp）也可以是chr start end（indel等）
        :return:
        """
        chr_data = defaultdict(list)
        with open(self.infile, 'r') as r:
            for line in r:
                if not re.match('#.*', line):
                    tmp = line.strip().split('\t')
                    chr_data[tmp[0]].append(int(tmp[1]))
        return chr_data

    def sliding_new(self):
        """
        换中方式去进行滑窗，之前的方法效率不高
        :return:
        """
        end = time.time()
        chr_lists = []
        win_data = {}
        step, i, x = 0, 0, 0
        s_step = self.step
        chr_data = self.read_vcf_file()
        lines = []
        for key in chr_data.keys():
            pos = chr_data[key]
            pos.sort()
            for jj in pos:
                lines.append([key, jj])
        end1 = time.time()
        print end1 - end
        for l in range(len(lines)):
            item = lines[l]
            if item[0] not in win_data.keys():
                chr_lists.append(item[0])
                win_data[item[0]] = []
                i = 0
                x = 0
                step = s_step
            if int(item[1]) <= step:
                x += 1
            else:
                n = int(item[1]) / s_step
                win_data[item[0]].append(x)
                i += 1
                if n > i:
                    m = i
                    for j in range(m, n):
                        step += s_step
                        x = 0
                        win_data[item[0]].append(x)
                        i += 1
                x = 1
                step += s_step
            if l < len(lines) - 1:
                tmp = lines[l + 1]
                if item[0] != tmp[0]:
                    win_data[item[0]].append(x)
        # noinspection PyBroadException
        try:
            win_data[item[0]].append(x)
        except:
            print "{}文件为空！".format(self.infile)
        end2 = time.time()
        # print win_data
        print end2 - end1
        outfile = os.path.join(self.output_dir, "win_{}_result.txt".format(self.step))
        with open(outfile, 'w') as w:
            w.write('#chrname\twin_result\n')
            for keys in self.sort_chr_sca(win_data):
                result = self.set_string_list(win_data[keys])
                w.write("{}\t{}\n".format(keys, ','.join(result)))

    def set_string_list(self, list__):
        string_list = []
        for n in list__:
            string_list.append(str(n))
        return string_list

    def sort_chr_sca(self, win_data):
        chr_list = []
        for key in win_data.keys():
            chr_list.append(key)
        try:
            chr_list.sort(key=lambda i: int(re.findall("\d+", i)[0]))
        except:
            chr_list.sort()
        return chr_list


if __name__ == '__main__':
    """
    用于对不同的位置信息文件进行划窗计算
    """
    parser = argparse.ArgumentParser(description="用于对不同的位置信息文件进行划窗计算。")
    parser.add_argument("-i", "--infile", type=str, help="输入含有位置信息的文件，必须保证前三列为chr,start,end！")
    parser.add_argument("-s", "--step", type=int, help="设置划窗的步长！", default=10000)
    parser.add_argument("-o", "--outfile", type=str, help="划窗结果文件所在路径")
    args = parser.parse_args()
    a = WindowSliding(args.infile, args.step, args.outfile)
    a.sliding_new()
