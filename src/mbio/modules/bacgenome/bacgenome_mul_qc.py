# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/8'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir,link_file
import pandas as pd
import time


class BacgenomeMulQcModule(Module):
    """
    微生物基因组多样本质控，二代数据: PE文库或MP文库；三代数据：pacbio数据及nanorepore数据的fastq
    """

    def __init__(self, work_id):
        super(BacgenomeMulQcModule, self).__init__(work_id)
        option = [
            {"name": "fastq_dir", "type": "infile", "format": "bacgenome.simple_dir", "required": True},
            {"name": "qc_a_list", "type": "infile", "format": "bacgenome.simple_file", "required": False},
            {"name": "qc_b_list", "type": "infile", "format": "bacgenome.simple_file", "required": False},
            {"name": "third_list", "type": "infile", "format": "bacgenome.simple_file", "required": False},
            {"name": "stat_list", "type": "infile", "format": "bacgenome.simple_file", "required": False},
            {"name": "phix_tool", "type": "string", "default": "bwa", "choose": ["bwa", "bowtie"]},
            {"name": "stat_new", "type": "outfile", "format": "bacgenome.simple_file", "required": False}
        ]
        self.add_option(option)
        self.qc_a = self.add_module("bacgenome.qc_a")  # fastp质控，包含原始数据统计，去phix，质控，质控后统计
        self.qc_b = self.add_module("bacgenome.qc_b")  # 旧的质控方法
        self.hiseq_stat = self.add_module("bacgenome.hiseq_seq_stat")  # 质控后的统计，对于不需质控的数据直接执行此步
        self.qc_third = self.add_module("bacgenome.pacbio_stat")  # 对三代数据进行长度统计
        self.run_modules = []  # 此list下的任务结束后set_output()
        self.sample_major_dir = os.path.join(self.work_dir, "major_qc")  # 同一样本多个文库下数据量最大的那组，用于三代数据校正和三代组装
        # self.sample_qc_tmp_dir = os.path.join(self.work_dir, "qc_tmp")  # 各子模块生成的结果整理至此
        self.sample_qc_result_dir = os.path.join(self.work_dir, "qc_result")  # 统一整理结果，分为[data, qc_stat, raw_stat, third_stat]四个目录
        self.stat_list_info = ""
        self.clean_data = ""
        self.lib_set = set()  # 一共有哪几种library
        self.all_samp,self.draft_samp,self.pe_samp,self.chr_samp, self.chr_pure, self.chr_mix = [[]] * 6

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("qc_a_list").is_set:
            if not self.option("qc_b_list").is_set:
                if not self.option("stat_list").is_set:
                    if not self.option("third_list").is_set:
                        raise OptionError("请提供至少一种list文件")
        if self.option("qc_a_list").is_set:
            self.run_modules.append(self.qc_a)
        if self.option("qc_b_list").is_set:
            self.run_modules.append(self.qc_b)
        if self.option("stat_list").is_set:
            self.run_modules.append(self.hiseq_stat)
        if self.option("third_list").is_set:
            self.run_modules.append(self.qc_third)
        self.store_samples_for_api()  # 保存哪些样本是扫描图，哪些样本是完成图，哪些样本有pe数据，哪些样本只有三代数据
        return True

    def store_samples_for_api(self):
        have_pe_set = set()
        have_pe_set |= self.store_samples(self.option("qc_a_list"))
        have_pe_set |= self.store_samples(self.option("qc_b_list"))
        have_pe_set |= self.store_samples(self.option("stat_list"))
        have_comp_set = self.store_samples(self.option("third_list"))
        self.pe_samp = list(have_pe_set)
        self.chr_samp = list(have_comp_set)
        self.all_samp = list(have_pe_set | have_comp_set)
        self.draft_samp = list(have_pe_set - have_comp_set)
        self.chr_pure = list(have_comp_set - have_pe_set)
        self.chr_mix = list(have_comp_set & have_pe_set)

    def store_samples(self, sample_list):
        if sample_list.is_set:
            data = pd.read_table(sample_list.prop["path"], header=None)
            sample_set = set(data[0].apply(str))
            self.lib_set |= set(data[1])
            return sample_set
        else:
            return set()

    def run(self):
        super(BacgenomeMulQcModule, self).run()
        self.on_rely(self.run_modules, self.set_output)
        if self.option("qc_a_list").is_set:
            self.get_list_info(self.option("qc_a_list").prop["path"])
            self.run_qc_a()
        if self.option("qc_b_list").is_set:
            self.get_list_info(self.option("qc_b_list").prop["path"])
            self.run_qc_b()
        if self.option("stat_list").is_set:
            self.get_list_info(self.option("stat_list").prop["path"])
            self.get_stat_list()
            self.run_hiseq_stat()
        if self.option("third_list").is_set:
            self.run_qc_third()

    def get_stat_list(self):
        self.stat_list = os.path.join(self.work_dir, "stat.txt")
        data = pd.read_table(self.option("stat_list").prop["path"], header=None)
        # data["sample_lib"] = data[0] + "_" + data[1] + "_" + data[2]
        data["sample_lib"] = data.apply(lambda x: x[6][:-5], axis=1)
        data[[6, "sample_lib", 5]].to_csv(self.stat_list, sep="\t", header=False, index=False)

    def get_list_info(self, list_file):
        data = pd.read_table(list_file, header=None)
        data["sample_lib"] = data.apply(lambda x: x[6][:-5], axis=1)
        data = data.set_index(["sample_lib"])
        if self.stat_list_info:
            self.stat_list_info = pd.concat([self.stat_list_info, data])
        else:
            self.stat_list_info = data


    def run_qc_a(self):
        self.qc_a.set_options({
            "list": self.option("qc_a_list"),
            "fastq_dir": self.option("fastq_dir"),
            "phix_tool": self.option("phix_tool")
        })
        self.qc_a.run()

    def run_qc_b(self):
        self.qc_b.set_options({
            "list": self.option("qc_b_list"),
            "fastq_dir": self.option("fastq_dir"),
            "phix_tool": self.option("phix_tool"),
        })
        self.qc_b.run()

    def run_hiseq_stat(self):
        self.hiseq_stat.set_options({
            "list": self.stat_list,
            "fastq_dir": self.option("fastq_dir")
        })
        self.hiseq_stat.run()

    def run_qc_third(self):
        self.qc_third.set_options({
            "third_list": self.option("third_list"),
            "fastq_dir": self.option("fastq_dir")
        })
        self.qc_third.run()

    def combine_stat(self):
        self.stat_info = pd.concat([self.qc_a.stat_info, self.qc_b.stat_info])

    def parse_stat(self, file_path):
        data = pd.read_table(file_path)
        self.stat_info = self.stat_info.append(data)

    def add_more_info(self, df):
        prefix = df["Sample_lib"]
        info = self.stat_list_info.loc[prefix]
        if isinstance(info, pd.DataFrame):
            info = info.iloc[0]
        df["Sample Name"] = info[0]
        df["Library"] = info[1]
        df["Library Name"] = info[2]
        df["Insert Size(bp)"] = info[3]
        df["Read Len(bp)"] = info[4]
        return df

    def write_stat(self):
        # 需要优化，写到work_dir里面的数据，和写在output里面的数据
        # 需要三种文件：1.提供下载的2.提供导表的3.提供流程用的
        if len(self.stat_info) == 0:
            return True
        stat_path = os.path.join(self.output_dir, "stat")
        if not os.path.isdir(stat_path):
            os.mkdir(stat_path)
        clean_data = self.stat_info[self.stat_info["type"] == "clean"].drop("type", axis=1)
        raw_data = self.stat_info[self.stat_info["type"] == "raw"].drop("type", axis=1)
        self.clean_data = clean_data.apply(self.add_more_info, axis=1)
        self.clean_data.reindex(columns=["Sample Name", "Library", "Library Name", "Insert Size(bp)", "Read Len(bp)",
                                    "pair reads(#)", "total bases(bp)", "Q20(%)", "Q30(%)"]).\
            to_csv(stat_path + "/clean_statistics.xls", sep="\t", header=True, index=False)
        self.clean_data.reindex(columns=["Sample Name", "Library", "Library Name", "Insert Size(bp)", "Read Len(bp)",
                                    "pair reads(#)", "total bases(bp)", "Q20(%)", "Q30(%)", "Sample_lib"]).\
            to_csv(stat_path + "/clean_statistics_mongo.xls", sep="\t", header=True, index=False)
        if len(raw_data) > 0:
            tmp = raw_data.apply(self.add_more_info, axis=1)
            tmp.reindex(columns=["Sample Name", "Library", "Library Name", "Insert Size(bp)", "Read Len(bp)",
                                        "pair reads(#)", "total bases(bp)", "Q20(%)", "Q30(%)"]).\
                to_csv(stat_path + "/raw_statistics.xls", sep="\t", header=True, index=False)
            tmp.reindex(columns=["Sample Name", "Library", "Library Name", "Insert Size(bp)", "Read Len(bp)",
                                        "pair reads(#)", "total bases(bp)", "Q20(%)", "Q30(%)", "Sample_lib"]).\
                to_csv(stat_path + "/raw_statistics_mongo.xls", sep="\t", header=True, index=False)

    def write_third_stat(self, path):
        data = pd.read_table(path)
        data_ref = pd.read_table(self.option("third_list").prop["path"], header=None)
        data["Library"] = data.apply(self.parse_third, axis=1, ref=data_ref)
        data.reindex(columns=["Sample", "Library", "Total Reads No.", "Total Bases (bp)", "Largest (bp)", "Average Len (bp)"]).\
            to_csv(path, sep="\t", header=True, index=False)

    def parse_third(self, df, ref):
        name = df[0]
        library = ref[ref[0]==name].iloc[0, 1]
        return library

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.combine_stat()
        fastx_path = os.path.join(self.output_dir, "fastx")
        data_path = os.path.join(self.output_dir, "cleandata")
        if not os.path.isdir(fastx_path):
            os.mkdir(fastx_path)
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        link_dir(self.qc_a.output_dir, self.output_dir)
        link_dir(self.qc_b.output_dir, self.output_dir)  # 会覆盖qc_a的q20q30统计结果，但是目前会重新生成这张表
        link_dir(self.hiseq_stat.output_dir, self.output_dir)
        if self.option("third_list").is_set:
            link_dir(self.qc_third.output_dir, self.output_dir + "/len")
            self.write_third_stat(self.output_dir+ "/len/statistics.xls")  # 重新改写statistics.xls，添加Library
        # 只stat的需要单独写
        for tool in self.hiseq_stat.stat_tools:
            files = os.listdir(tool.output_dir)
            for file in files:
                file_path = os.path.join(tool.output_dir, file)
                if file == "qc_stat.xls":
                    self.parse_stat(file_path)
                else:
                    file_prefix = file[:-15]
                    if file[-14] == "1":
                        new_file = file_prefix + "_l.raw_fastxstat"
                    else:
                        new_file = file_prefix + "_r.raw_fastxstat"
                    new_file_path = os.path.join(fastx_path, new_file)
                    link_file(file_path, new_file_path)
        self.write_stat()
        self.logger.info("设置qc结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BacgenomeMulQcModule, self).end()
