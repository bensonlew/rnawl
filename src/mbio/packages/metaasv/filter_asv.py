# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

"""
功能：对输入的asv表 进行过滤
"""
import re
import json
import argparse
from optparse import OptionParser
from mbio.packages.metaasv.common_function import filter_asv
from biocluster.config import Config

class FilterASV(object):
    """
    {"species_filter" : [{"type":"0","level_id"：1, "value":"d__Bacteria"}],
    "sample_filter" : {"type":"0", "sample_num":"3", "reads_num":"0.05"},
    "reads_filter" : {"type":"0","reads_num":"0.05"},
    "func_filter" : "0",
    "set_list" : [{"type":"0","set_id": "dasdfaf"}}
    ]
    type类型：1代表去除，0代表保留
    """
    def __init__(self, otu_table, filter_json,db_version):
        super(FilterASV, self).__init__()
        self.otu_table = otu_table ##丰度表为相对丰度表
        self.otu_sample_dict = self.extract_info(self.otu_table)
        self.filter_json = filter_json
        self.db_version = db_version
        with open(otu_table, 'rb') as r:
            header = r.readline()
            self.otu_json = r.readlines()
        self.LEVEL = {
            1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
            6: "f__", 7: "g__", 8: "s__", 9: "asv"
        }

    def filter_asv_table(self):
        """
        筛选 ASV表中的物种，样本和reads
        params:
        """
        # print(self.filter_json)
        filter = json.loads(self.filter_json)
        # print(filter)
        filter_json = self.byteify(filter)
        # print(filter_json)
        sp_condition_keep = []
        sp_condition_remove = []
        sam_condition = []
        reads_condition = []
        func_filter = []
        set_list = []
        if "species_filter" in filter_json.keys():
            species_list = filter_json["species_filter"]
            for condition in species_list:
                condition_level = condition["value"]
                condition_type = condition["type"]
                condition['value'] = self.escapeExprSpecialWord(condition_level)
                if condition_type in ["0", 0]:
                    sp_condition_keep.append(condition)
                if condition_type in ["1", 0]:
                    sp_condition_remove.append(condition)
        if "sample_filter" in filter_json.keys():
            sam_condition.append(filter_json["sample_filter"])
        if "reads_filter" in filter_json.keys():
            reads_condition.append(filter_json["reads_filter"])
        if "func_filter" in filter_json.keys():
            func_filter.append(filter_json["func_filter"])
        if "asv_filter" in filter_json.keys():
            set_list= filter_json["asv_filter"]
        if sp_condition_keep:
            # print(sp_condition_keep)
            self.otu_json = self.keep_species(sp_condition_keep)
        if sp_condition_remove:
            self.otu_json = self.remove_species(sp_condition_remove)
        if sam_condition:
            for i in sam_condition: #以后可能会出现多个条件的情况
                print(i)
                self.filter_samples(i)
        if reads_condition:
            for i in reads_condition:#以后可能会出现多个条件的情况
                print(i)
                self.filter_reads(i)
        if func_filter:
            for i in func_filter:#以后可能会出现多个条件的情况
                print(i)
                type = i['type']
                self.otu_json = self.filter_function(type)
        if set_list:
            for i in set_list:
                self.otu_json = self.filter_asv_json(i)
        return self.otu_json

    def extract_info(self, otu_table):
        """
        解析ASV表，将ASV表解析成一个二维字典 相对丰度表
        例如 dict[ASV1][sample1] = 0.025
        代表sample1有20条序列在ASV1里
        """
        info_dict = {}
        with open(otu_table, 'r') as r:
            line = r.next().strip().split("\t")
            sample_name = line[1:]
            for line in r:
                sample_dict = {}
                line = line.strip().split("\t")
                otu_name = line[0]
                for i in range(len(sample_name)):
                    sample_dict[sample_name[i]] = float(line[i+1])
                info_dict[otu_name] = sample_dict
        return info_dict

    def escapeExprSpecialWord(self, words):
        """
        转义正则表达式中的特殊字符
        :params words: 需要转义的字符串
        :return : 转义之后的字符串
        """
        fbsArr = ["\\", "$", "(", ")", "*", "+", ".", "[", "]", "?", "^", "{", "}", "|"]
        for key in fbsArr:
            if words.find(key) >= 0:
                words = words.replace(key, "\\" + key)
        return words

    def keep_species(self, conditions, fuzzy=False):
        """
        保留特定物种的ASV, 多条件时选交集
        :params conditions: 条件列表
        :params fuzzy: 是否模糊匹配(不区分大小写)
        :return : 筛选出来的ASV列表
        """
        temp_otus = []
        for cond in conditions:
            if not fuzzy:
                cond['pattern'] = re.compile('^{}$'.format(cond['value'].strip()))
            else:
                cond['pattern'] = re.compile('{}'.format(cond['value'].strip()), flags=re.IGNORECASE)
        for i in self.otu_json:
            for cond in conditions:
                level = int(cond["level_id"]) - 1
                sp_name = re.split(r';', (re.split(r'\t', i, maxsplit=1)[0]))
                sp_name = sp_name[level].strip()
                # print(sp_name)
                if cond['pattern'].search(sp_name):
                    temp_otus.append(i)
                    break
        return temp_otus

    def remove_species(self, conditions, fuzzy=False):
        """
        去除特定物种的ASV, 多条件时，任何条件去除即去除
        :params conditions: 条件列表
        :params fuzzy: 是否模糊匹配(不区分大小写)
        :return : 筛选出来的ASV列表
        """
        temp_otus = []
        for cond in conditions:
            if not fuzzy:
                cond['pattern'] = re.compile('^{}$'.format(cond['value'].strip()))
            else:
                cond['pattern'] = re.compile('{}'.format(cond['value'].strip()), flags=re.IGNORECASE)
        for i in self.otu_json:
            for cond in conditions:
                level = int(cond["level_id"]) - 1
                sp_name = re.split(r';', (re.split(r'\t', i, maxsplit=1)[0]))
                sp_name = sp_name[level].strip()
                if cond['pattern'].search(sp_name):
                    break
            else:
                temp_otus.append(i)
        return temp_otus

    def filter_samples(self, my_json):
        """
        保留或者去除至少在x个样本中相对丰度都小于y的物种(ASVs)
        """
        tmp_list = self.otu_json[:]
        my_c = 0
        for line in self.otu_json:
            otu = line.split("\t")[0]  # asv即是第一列，也是self.otu_sample_dict的第一个key值
            count = 0
            type = my_json['type']
            reads_number = my_json["reads_num"]
            for sp in self.otu_sample_dict[otu]:
                if sp in ["Percent"]:
                    continue
                else:
                    if float(self.otu_sample_dict[otu][sp]) < float(reads_number):
                        count += 1 ###  计算样本个数
            if count >= int(my_json["sample_num"]):
                if type in ["0", 0]:
                    continue
                elif type in ["1", 1]:
                    tmp_list.remove(line)
            else:
                if type in ["0", 0]:
                    tmp_list.remove(line)
                elif type in ["1", 1]:
                    continue
                my_c += 1
        print(len(tmp_list))
        self.otu_json = tmp_list[:]
        print(len(self.otu_json))

    def filter_reads(self, my_json):
        """
        保留或者去除序列数总和小于x的物种(ASV)
        """
        tmp_list = self.otu_json[:]
        for line in self.otu_json:
            otu = line.strip().split("\t")[0]  # asv即是第一列，也是self.otu_sample_dict的第一个key值
            summary = 0
            type = my_json['type']
            reads_num = my_json["reads_num"]
            for sp in self.otu_sample_dict[otu]: ##遍历样本
                if sp in ["Percent"]:
                    summary = float(self.otu_sample_dict[otu][sp])
            # print("summary: {}".format(summary))
            if summary <= float(reads_num):
                if type in ["0", 0]:
                    pass
                elif type in ["1", 1]:
                    tmp_list.remove(line)
            else:
                if type in ["0", 0]:
                    tmp_list.remove(line)
                elif type in ["1", 1]:
                    pass
        print(len(tmp_list))
        self.otu_json = tmp_list[:]
        print(len(self.otu_json))

    def filter_function(self, type):
        """
        保留或者去除 比对到叶绿体（Chloroplast）和线粒体（Mitochondrial）序列的ASVs
        叶绿体：o__Chloroplast, 线粒体：f__Mitochondria
        """
        temp_otus = self.otu_json[:]
        tmp_list = []
        remove_list = ["o__Chloroplast", "f__Mitochondria"]
        type = str(type)
        for i in self.otu_json:
            # sp_name = re.split(r';', (re.split(r'\t', i, maxsplit=1)[0]))
            sp_name = re.split(r'\t', i, maxsplit=1)[0]
            if type in ["0", 0]:
                for cond in remove_list:
                    if cond in sp_name:
                        tmp_list.append(i)
            elif type in ["1", 1]:
                for cond in remove_list:
                    if cond in sp_name:
                        temp_otus.remove(i)
        print(len(temp_otus))
        if type in ["1", 1]:
            return temp_otus
        else:
            return tmp_list

    def filter_asv_json(self,my_json):
        """
        过滤asv集
        :return:
        """
        # print("===+++++:"+ my_json)
        input_json = my_json
        set_id = input_json["set_id"]
        asv_list = filter_asv(set_id,self.db_version)
        type = input_json["type"]
        temp_otus = []
        type = str(type)
        for i in self.otu_json:
            sp_name = re.split(r';', (re.split(r'\t', i, maxsplit=1)[0]))[-1]##取出来asv
            if type in ["0", 0]:
                if sp_name in asv_list:
                    temp_otus.append(i)
            elif type in ["1", 1]:
                if sp_name in asv_list:
                    continue
                else:
                    temp_otus.append(i)
        return temp_otus

    def byteify(self, input):
        if isinstance(input, dict):
            return {self.byteify(key): self.byteify(value) for key, value in input.iteritems()}
        elif isinstance(input, list):
            return [self.byteify(element) for element in input]
        elif isinstance(input, unicode):
            return input.encode("utf-8")
        else:
            return input

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='input', type=str, help='input dir file')
    parser.add_argument('-j', metavar='filter_json', type=str, help='filter_json')
    parser.add_argument('-o', metavar='output', type=str, help='output asv table')
    parser.add_argument('-d', metavar='db', type=str, help='mongodb version')
    args = parser.parse_args()
    input_file = args.i
    filter_json = args.j
    output = args.o
    if args.d:
        Config().DBVersion = args.d
    filter = FilterASV(input_file, filter_json,args.d)
    asv_table_list = filter.filter_asv_table()
    with open(input_file, 'r') as f, open(output, "w") as w:
        header = f.readline()
        w.write(header)
        for asv in asv_table_list:
            w.write(asv)