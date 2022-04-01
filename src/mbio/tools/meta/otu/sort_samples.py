# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import os
import shutil
import linecache
import numpy as np
import pandas as pd
from collections import defaultdict
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from multiprocessing import Process
from mbio.packages.taxon.mask_taxon import mask_taxon  #add guanqing 20180807
import re

class SortSamplesAgent(Agent):
    """
    传入一个group表，以及是否进行样本合并的参数生成一张OTU表并对并依照group表OTU表进行筛选合并
    """
    def __init__(self, parent):
        super(SortSamplesAgent, self).__init__(parent)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU文件
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 输入的group表
            {"name": "method", "type": "string", "default": ""},  # 样本的合并方式, ""为不进行合并
            {"name": "out_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的结果OTU表
            {"name": "level_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的结果OTU表(百分比）
            {"name": "others", "type": "float", "default": ""},  # 组成分析中用于将丰度小于0.01/其它的物种归为others
            {"name": "top_n", "type": "int", "default": 0}  # 获取总丰度前 n 的物种
        ]
        self.add_option(options)
        self.step.add_steps("sort_samples")
        self.on('start', self.start_sort_samples)
        self.on('end', self.end_sort_samples)

    def start_sort_samples(self):
        self.step.sort_samples.start()
        self.step.update()

    def end_sort_samples(self):
        self.step.sort_samples.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option("in_otu_table").is_set:
            raise OptionError("输入的OTU文件不能为空", code="32705301")
        if self.option("method"):
            if self.option("method") not in ["", "no", "none", "No", "None", None, "average", "sum", "middle"]: # add middle by zhouxuan 20161205
                raise OptionError("参数method设置错误！", code="32705302")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["taxa.table.xls", "xls", "各样本物种丰度结果表"], # modified by hongdongxuan 20170323 otu_otu 改为taxa.table
            ["taxa.percents.table.xls", "xls", "各样本物种相对丰度结果表"]  # add by wangzhaoyue 2017.03.06 modified by hongdongxuan 20170323 level_otu_table 改为taxa.precents.table
        ])
        super(SortSamplesAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "10G"  # modified by guhaiodng 20171101


class SortSamplesTool(Tool):
    def __init__(self, config):
        super(SortSamplesTool, self).__init__(config)
        self.logger.info("SortSamples读取分组表开始")
        samples = list()
        with open(self.option("group_table").prop["path"], "rb") as r:
            line = r.next()
            for line in r:
                line = line.rstrip().split("\t")
                samples.append(line[0])
        self.samples = samples
        self.logger.info("SortSamples读取分组表结束")

    def filter_samples(self):
        no_zero_otu = os.path.join(self.work_dir, "otu.nozero")

        #self.option("in_otu_table", self.work_dir + "/tmp_mask_otu.xls")  #guanqing
        self.option("in_otu_table").sub_otu_sample(self.samples, no_zero_otu)
        no_zero_otu_new = self.work_dir + "/tmp_mask_otu.xls"
        self.name_to_name = mask_taxon(no_zero_otu, no_zero_otu_new) #guanqing 20180807

        #return no_zero_otu
        return no_zero_otu_new

    ##add  guanqing 20180807
    def dashrepl(self, matchobj):
        """
        add func by guhaidong 20171025
        """
        return self.name_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result):
        """
        add func by guhaidong 20171025
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                # line = i.strip()
                new_line = re.sub(r"(name\d+)", self.dashrepl, i)
                w.write(new_line)

    def cat_samples(self, otu, method):
        sample_group = dict()  # 一个样本是属于哪个group的
        index_sample = dict()  # 一个OTU表中第几列属于哪个样本
        group_sample_num = defaultdict(int)  # 一个分组里面有多少的样本
        cat_otu_path = os.path.join(self.work_dir, "cat_otu.xls")
        with open(self.option("group_table").prop["path"], "rb") as r:
            line = r.next()
            for line in r:
                line = line.rstrip().split("\t")
                sample_group[line[0]] = line[1]
                group_sample_num[line[1]] += 1
        with open(otu, "rb") as r, open(cat_otu_path, 'wb') as w:
            group_list = list()
            for v in sample_group.values():
                group_list.append(v)
                group_list = list(set(group_list))
                # l = len(group_list) #zx

            line = r.next().rstrip().split("\t")
            for i in range(len(line)):
                index_sample[i] = line[i]

            w.write(index_sample[0] + "\t")
            w.write("\t".join(group_list) + "\n")
            for line in r:
                line = line.rstrip().split("\t")
                num = defaultdict(int)
                # middle_num = defaultdict(int)
                tmp = list()
                list1 = []  # add 2 lines by zhouxuan 20161205
                mid_num = dict()

                w.write(line[0] + "\t")
                for i in range(1, len(line)):
                    num[sample_group[index_sample[i]]] += float(line[i])
                for m in group_list:  # add 12 lines by zhouxuan 20161205
                    for i in range(1, len(line)):
                        if sample_group[index_sample[i]] == m:
                            list1.append(float(line[i]))
                            if len(list1) == group_sample_num[m]:
                                list1.sort()
                                yu = int(group_sample_num[m]) % 2
                                index = int(int(group_sample_num[m]) / 2)
                                if yu == 0:
                                    mid_num[m] = int(round((int(list1[index-1]) + int(list1[index]))/2))
                                    list1 = []
                                else:
                                    mid_num[m] = list1[index]
                                    list1 = []

                if method == "sum":
                    for g in group_list:
                        tmp.append(str(num[g]))
                if method == "average":
                    for g in group_list:
                        avg = int(round(num[g] / group_sample_num[g]))
                        tmp.append(str(avg))
                if method == "middle":  # add 3 line by zhouxuan 20161205
                    for g in group_list:
                        tmp.append(str(mid_num[g]))
                """
                else :
                    if method == "average":
                        for g in group_list:
                            avg = int(round(num[g] / group_sample_num[g]))
                            tmp.append(str(avg))
                    else method == "middle": # add 6 line by zhouxuan 20161205
                        for g in group_list:
                            list1.sort()
                            index = int(round(len(line)-1 / 2))
                            middle_num = list1[index]
                            tmp.append(str(middle_num))
                """

                w.write("\t".join(tmp))
                w.write("\n")
        return cat_otu_path

    def percent(self, origin_file):   # add by wangzhaoyue 2017.03.06
        tmp = np.loadtxt(origin_file, dtype=np.str, delimiter="\t")
        data = tmp[1:, 1:].astype(np.float)
        data[np.isnan(data)] = 0.0 ## 修改otu表中某些样本为空的情况
        array_data = data / data.sum(0)
        array = np.c_[tmp[1:, 0], array_data]
        array_all = np.row_stack((tmp[0, :], array))
        s = np.char.encode(array_all, 'utf-8')
        #percent_file = os.path.join(self.output_dir, "taxa.percents.table.xls")
        percent_file = os.path.join(self.work_dir, "tmp_taxa.percents.table.xls")
        np.savetxt(percent_file, s, fmt='%s',  delimiter='\t', newline='\n')
        return percent_file

    def get_others(self, percent_file):  # add by zhujuan 2017.10.12
        df = pd.DataFrame(pd.read_table(percent_file, sep='\t', index_col=0))
        new_df = df.ix[list((df > self.option("others")).any(axis=1))]
        new_df2 = new_df.copy()
        others = df.ix[list((df < self.option("others")).all(axis=1))]
        if len(others) > 0 :
            new_df2.loc["others"] = others.apply(lambda x: x.sum(), axis=0)
        #other = os.path.join(self.output_dir, "taxa.percents.table.xls")
        other = os.path.join(self.work_dir, "tmp_taxa.percents.table.xls")
        new_df2.to_csv(other, sep="\t", encoding="utf-8")

    def get_top(self, infile, top_n=[]):
        '''
        获取总丰度前 n 的物种
        如果没指定top_n的list，就先计算
        '''
        tb = pd.read_csv(infile, sep='\t', index_col=0)
        index_name = tb.index.name
        if not any(top_n):
            row_sum = tb.sum(1).sort_values(ascending=False)
            top_n = row_sum.index[:self.option('top_n')]
        tb = tb.loc[top_n]
        tb.index.name = index_name
        tb.to_csv(infile, sep='\t')
        return infile, top_n

    def run(self):
        super(SortSamplesTool, self).run()
        final_otu = self.filter_samples()
        if self.option("method") in ["average", "sum", "middle"]:
            final_otu = self.cat_samples(final_otu, self.option("method"))
        #out_otu = os.path.join(self.output_dir, "taxa.table.xls")  # modified by hongdongxuan 20170323 otu_otu 改为taxa.table
        out_otu_new = os.path.join(self.output_dir, "taxa.table.xls")
        #shutil.copy2(final_otu, out_otu)
        self.add_taxon(final_otu, out_otu_new) #guanqing
        #self.percent(out_otu)
        #self.percent(final_otu)  #guanqing
        #final_level_percents = self.percent(out_otu)
        final_level_percents = self.percent(final_otu)
        if self.option("others") != "":
            self.get_others(final_level_percents)

        final_level_percents_new = os.path.join(self.output_dir, "taxa.percents.table.xls")  #guanqing
        self.add_taxon(final_level_percents, final_level_percents_new)  #guanqing
        #self.option("level_otu_table").set_path(final_level_percents)   # add 3 lines by wangzhaoyue 2017.03.06
        if self.option('top_n'):
            out_otu_new, top_n = self.get_top(out_otu_new)
            final_level_percents_new, _ = self.get_top(final_level_percents_new, top_n)

        self.option("level_otu_table").set_path(final_level_percents_new)
        #self.option("out_otu_table").set_path(out_otu)
        self.option("out_otu_table").set_path(out_otu_new)
        self.end()
