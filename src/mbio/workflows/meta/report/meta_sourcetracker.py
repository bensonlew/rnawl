# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan' 20170503

"""微生物组成来源比例分析模块"""
import os
import json
import shutil
import datetime
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import re
import pandas as pd
from pandas import Series,DataFrame
from mainapp.models.mongo.public.meta.meta import Meta


class MetaSourcetrackerWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetaSourcetrackerWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "level", "type": "string", "default": "9"},  # 输入的OTU level
            {"name": "map_detail", "type": "infile", "format": "meta.otu.group_table"},  # 输入的map_detail 示例如下(map文件后续导表)
            {"name": "group_detail", "type": "string"},
            {"name": "second_group_detail", "type": "string"},
            {"name": "add_Algorithm", "type": "string", "default": ""},
            {"name": "s", "type": "string", "default": "1"},  # OTU筛选参数
            {"name": "meta_sourcetracker_id", "type": "string"}, #主表的id
            {"name": "update_info", "type": "string"}
            # {"name": "source", "type": "string"}
            # {"A":["578da2fba4e1af34596b04ce","578da2fba4e1af34596b04cf","578da2fba4e1af34596b04d0"],"B":["578da2fba4e1af34596b04d1","578da2fba4e1af34596b04d3","578da2fba4e1af34596b04d5"],"C":["578da2fba4e1af34596b04d2","578da2fba4e1af34596b04d4","578da2fba4e1af34596b04d6"]}
            # {"name": "method", "type": "string", "default": ""}  # 聚类方式， ""为不进行聚类
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.meta_sourcetracker = self.add_tool("meta.beta_diversity.meta_sourcetracker")
        self.sort_all_samples = self.add_tool("meta.otu.sort_samples")
        self.sort_source_samples = self.add_tool("meta.otu.sort_samples")
        self.sort_sink_samples = self.add_tool("meta.otu.sort_samples")
        self.all_sort = ''
        self.source_group_file = ''
        self.sink_group_file = ''
        self.math_sample_name = ''
        self.sample_name = dict()

    def check_options(self):  # 2016.12.1 zhouxuan
        if self.option('add_Algorithm') not in ['average', 'middle', 'sum', ""]:
            raise OptionError('错误的层级聚类方式：%s', variables=(self.option('add_Algorithm')), code="12702401")

    def judge(self):
        """
        判断map文件是否正确，同时生成后续要用的两个小的分组文件
        :return:
        """
        self.source_group_file = os.path.join(self.work_dir, "source_group")
        self.sink_group_file = os.path.join(self.work_dir, "sink_group")
        samples = []
        self.sink_label = []
        self.sample_name = dict()
        with open(self.option('map_detail').prop['path'], "rb") as r, open(self.source_group_file, 'a') as e, open(self.sink_group_file, 'a') as k:
            line = r.next()
            k.write("#Sample\tsink\n")
            e.write("#Sample\tsource\n")
            for line in r:
                line = line.rstrip().split("\t")
                if re.match('[0-9].+', line[0]):
                    self.math_sample_name = "true"
                if line[0] not in samples:
                    samples.append(line[0])
                else:
                    raise OptionError('sink组和source组中不能存在同一个样本 %s', variables=(line[0]), code="12702402")
                if line[2] == 'sink':
                    k.write(line[0] + "\t" + line[1] + "\n")
                    self.sample_name[line[0]] = "N_" + line[0]
                    if line[1] not in self.sink_label:
                        self.sink_label.append(line[1])
                else:
                    e.write(line[0] + "\t" + line[1] + "\n")
                    self.sample_name[line[0]] = "N_" + line[0]
        self.logger.info(self.sample_name)
        if len(self.sink_label) == 1 or self.option('add_Algorithm') == '':
            self.run_sort_all_samples()
        else:
            self.run_sort_source_samples()

    def run_sort_all_samples(self):
        self.all_sort = "true"
        self.sort_all_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.option("map_detail")
        })
        self.sort_all_samples.on('end', self.run_meta_sourcetracker)
        self.sort_all_samples.run()

    def run_sort_source_samples(self):
        self.sort_source_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.source_group_file
        })
        self.sort_source_samples.on('end', self.run_sort_sink_samples)
        self.sort_source_samples.run()

    def run_sort_sink_samples(self):
        self.sort_sink_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.sink_group_file,
            "method": self.option('add_Algorithm')
        })
        self.sort_sink_samples.on('end', self.run_meta_sourcetracker)
        self.sort_sink_samples.run()

    def run_meta_sourcetracker(self):
        if self.all_sort == "true":
            qiime_table_path = self.reset_input_otu_table(self.sort_all_samples.option("out_otu_table").prop['path'])
            map_file_path = self.option('map_detail').prop['path']
        else:
            otu_table = self.add_otu_file(self.sort_source_samples.option("out_otu_table").prop['path'], self.sort_sink_samples.option("out_otu_table").prop['path'])
            map_file_path = self.get_new_map()
            qiime_table_path = self.reset_input_otu_table(otu_table)
        if self.math_sample_name == "true":
            fin_map = self.change_name(map_file_path)
        else:
            fin_map = map_file_path
        self.meta_sourcetracker.set_options({
            "otu_table": qiime_table_path,
            "map_table": fin_map,
            "s": self.option("s")
        })
        self.meta_sourcetracker.on('end', self.set_db)
        self.meta_sourcetracker.run()

    def reset_input_otu_table(self, old_otu_table_path):
        self.logger.info("开始重置OTU表")
        new_qiime_otu_table = os.path.join(self.work_dir, "finally_table.txt")
        with open(new_qiime_otu_table, "a") as n:
            firstline = "# QIIME-formatted OTU table" + "\n"
            n.write(firstline)
        if os.path.exists(old_otu_table_path):
            a = open(old_otu_table_path, "r")
            content = a.readlines()
            m = ''
            for f in content:
                if f.startswith("OTU ID") == True:
                    if self.math_sample_name != "true":
                        line = "#" + f
                    else:
                        f = f.rstrip("\n").split("\t")
                        for c in f[1:]:
                            if c in self.sample_name.keys():
                                m = m + "\t" + self.sample_name[c]
                            else:
                                m = m + "\t" + "N_" + c
                                self.sample_name[c] = "N_" + c
                        line = "#OTU ID" + m + "\n"
                        print line
                else:
                    line = f
                with open(new_qiime_otu_table, "a") as n:
                    n.write(line)
            a.close()
        else:
            self.logger.error("qiime的前一张OTU表不存在")
            self.set_error("qiime的前一张OTU表不存在", code="12702403")
        self.logger.info("重置OTU表结束")
        return new_qiime_otu_table

    def add_otu_file(self, source_file, sink_file):
        self.logger.info("开始合并source和sink")
        new_otu_file = os.path.join(self.work_dir, "add_table.txt")
        source = pd.read_table(source_file, header=0, sep="\t")
        sink = pd.read_table(sink_file, header=0, sep="\t")
        source_spe = source["OTU ID"]
        source_new_index = [i for i in source_spe]
        source.index = source_new_index
        e = set(source_new_index)
        sink_spe = sink["OTU ID"]
        sink_new_index = [i for i in sink_spe]
        sink.index = sink_new_index
        k = set(sink_new_index)
        spe = e & k
        fin_spe = [i for i in spe]
        source.reindex(index=fin_spe)
        sink.reindex(index=fin_spe)
        fin = pd.merge(source, sink)
        fin.to_csv(new_otu_file, sep="\t", index=False)
        self.logger.info("合并source和sink结束")
        return new_otu_file

    def change_name(self, file1):
        self.logger.info("self.math_sample_name {}".format(self.math_sample_name))
        self.logger.info("修改分组文件的样本名称")
        new_file = os.path.join(self.work_dir, "finally_map_table.txt")
        with open(file1, "rb") as r, open(new_file, 'a') as e:
            for line in r:
                if line.startswith("#") == True:
                    e.write(line)
                else:
                    line = line.rstrip("\n").split("\t")
                    e.write(str(self.sample_name[line[0]]) + "\t" + str(("\t").join(line[1:]) + "\n"))
        return new_file

    def revert_name(self, file1, file2):
        self.logger.info("修改结果文件中的样本名称")
        key_list = []
        value_list = []
        for key, value in self.sample_name.items():
            key_list.append(key)
            value_list.append(value)
        with open(file1, "rb") as r, open(file2, 'a') as e:
            for line in r:
                if line.startswith("S") == True:
                    e.write(line)
                else:
                    line = line.rstrip("\n").split("\t")
                    e.write(key_list[value_list.index(line[0])] + "\t" + str(("\t").join(line[1:])) + "\n")
        return file2

    def get_new_map(self):
        """
        sink分组进行分组求和时产生新的map表
        :return:
        """
        self.logger.info("sink分组进行分组求和时产生新的map表")
        new_map = os.path.join(self.work_dir, "group_map.txt")
        with open(self.option('map_detail').prop['path'], "rb") as r, open(new_map, 'a') as e:
            for line in r:
                line = line.rstrip("\n").split("\t")
                self.logger.info(len(line))
                if line[2] != "sink":
                    e.write(str(("\t").join(line)) + "\n")
            for c in self.sink_label:
                e.write(c + "\t" + "s_" + c + "\t" + "sink\n")  # 为了保证样本名称和分组名称的不一致性 所以强加s_
        self.logger.info("sink分组进行分组求和时产生新的map表的生成结束")
        return new_map

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        if self.math_sample_name != "true":
            file_path = self.meta_sourcetracker.output_dir + "/sink_predictions.txt"
            stdev_file_path = self.meta_sourcetracker.output_dir + "/sink_predictions_stdev.txt"
        else:
            file_path = self.revert_name(self.meta_sourcetracker.output_dir + "/sink_predictions.txt",
                                         os.path.join(self.work_dir, "sink_predictions.txt"))
            stdev_file_path = self.revert_name(self.meta_sourcetracker.output_dir + "/sink_predictions_stdev.txt",
                                               os.path.join(self.work_dir, "sink_predictions_stdev.txt"))
        api_otu = self.api.meta_sourcetracker
        api_otu.add_sg_sourcetracker_detail(self.option("meta_sourcetracker_id"), file_path=file_path,
                                            stdev_file_path=stdev_file_path,
                                            name_1="sink_predictions.txt",
                                            name_2="sink_predictions_stdev.txt")
        self.end()

    def end(self):
        if self.math_sample_name != "true":
            try:
                os.link(self.meta_sourcetracker.output_dir + "/sink_predictions.txt",
                        self.output_dir + "/sink_predictions.txt")
                os.link(self.meta_sourcetracker.output_dir + "/sink_predictions_stdev.txt",
                        self.output_dir + "/sink_predictions_stdev.txt")
            except Exception as e:
                self.logger.error("copying sink_predictions.txt failed{}".format(e))
                self.set_error("copying sink_predictions.txt failed", code="12702401")
                raise Exception("copying sink_predictions.txt failed{}".format(e))
        else:
            try:
                os.link(self.work_dir + "/sink_predictions.txt",
                        self.output_dir + "/sink_predictions.txt")
                os.link(self.work_dir + "/sink_predictions_stdev.txt",
                        self.output_dir + "/sink_predictions_stdev.txt")
            except Exception as e:
                self.logger.info("copying sink_predictions.txt failed{}".format(e))
                self.set_error("copying sink_predictions.txt failed", code="12702402")
                raise Exception("copying sink_predictions.txt failed{}".format(e))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "微生物组成来源比例分析"],
            ["./sink_predictions.txt", "txt", "相对贡献度表"],
            ["./sink_predictions_stdev.txt", "txt", "相对贡献度标准差表"]
        ])
        super(MetaSourcetrackerWorkflow, self).end()

    def run(self):
        self.judge()
        super(MetaSourcetrackerWorkflow, self).run()
