# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file
import pandas as pd
import os


class Qiime2DenoiseModule(Module):
    """
    qiime2 降噪模块
    """
    def __init__(self, work_id):
        super(Qiime2DenoiseModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},##fastq_dir文件夹
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列,自定义模式上传的fasta文件
            {'name': 'denoise_method', 'type': 'string'},  # 降噪方法选择
            {'name': 'database', 'type': 'string'},  # 数据库选择  用于选择参考数据库
            {"name": "truc_len", "type": "int", "default": 0},## 过滤长度的值
            {"name": "sample_type", "type": "string", "default": "multiple"},# 单样本还是多样本
            {"name": "max_ee", "type": "int", "default": 2},
            {"name": "trunc_q", "type": "int", "default": 0}, ##判断打断序列的大小，如果为0则不打断，即不质控
            {"name": "min_size", "type": "int", "default": 1},
            {"name": "min_reads", "type": "int", "default": 1},
        ]
        self.add_option(options)
        self.file_list = []

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("请输入fastq文件夹！")
        if self.option("denoise_method") in ['DADA2, Deblur']:
            raise OptionError("请输入降噪的方法！")

    def run_data2(self):
        """
        运行软件qiime2 data2_denoise
        :return:
        """
        if self.option("sample_type") in ['multiple']:
            self.data2 = self.add_module("metaasv.dada2_denoise")
        else:
            self.data2 = self.add_module("metaasv.data2_denoise")
        options = {
            "fastq_dir":self.option("fastq_dir"),
            "max_ee": self.option("max_ee"),
            "trunc_q": self.option("trunc_q"),
            }
        self.data2.set_options(options)
        self.data2.on("end",self.set_output)
        self.data2.run()

    def run_deblur(self):
        """
        运行软件qiime2 deblur_denoise
        :return:
        """
        if self.option("sample_type") in ['multiple']:
            self.deblur = self.add_module("metaasv.deblur_denoise_split")
        else:
            self.deblur = self.add_module("metaasv.deblur_denoise")
        options = {
            "fastq_dir":self.option("fastq_dir"),
            "database": self.option("database"),
            # "trim_length": self.option("trim_length"),
            }
        if self.option("ref_fasta").is_set:
            options["ref_fasta"] = self.option("ref_fasta")
        if self.option("truc_len") != 0:
            options["truc_len"] = self.option("truc_len")
        if self.option("min_reads") != 0:
            options["min_reads"] = self.option("min_reads")
        if self.option("min_size"):
            options["min_size"] = self.option("min_size")
        self.deblur.set_options(options)
        self.deblur.on("end",self.set_output)
        self.deblur.run()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        if self.option("denoise_method") in ["DADA2"]:
            link_dir(self.data2.output_dir, self.output_dir)
        else:
            link_dir(self.deblur.output_dir, self.output_dir)
        self.create_mongo_data()
        self.logger.info("设置结果文件目录成功！")
        if os.path.exists(os.path.join(self.output_dir, "stats.tsv")):
            os.remove(os.path.join(self.output_dir, "stats.tsv"))
        if os.path.exists(os.path.join(self.output_dir, "feature-table.biom")):
            os.remove(os.path.join(self.output_dir, "feature-table.biom"))
        self.end()

    def create_mongo_data(self):
        """
        将ASV_table.xls和降噪后的序列文件统计信息进行合并
        :return:
        """
        asv_table = pd.read_table(self.output_dir + "/ASV_table.xls", sep="\t", header=0)
        mongo_denoise_stat = self.output_dir + "/denoise_stat.xls"
        sample_list = list(asv_table.columns)[1:]
        sample_number_dict = {}
        sequence_number_dict = {}
        for sample in sample_list:
            asv_number = len(list(asv_table[sample][asv_table[sample] != 0.0]))
            if sample not in sample_number_dict:
                sample_number_dict[sample] = asv_number
            sequence_number = asv_table[sample][asv_table[sample] != 0.0].values.sum()
            if sample not in sequence_number_dict:
                sequence_number_dict[sample] = sequence_number
        if self.option("denoise_method") in ["DADA2"]:
            denoise_path = self.output_dir + "/DADA2_sequence_info.txt"
            with open(denoise_path, 'r') as f, open(mongo_denoise_stat, 'w') as w:
                lines = f.readlines()
                w.write("Sample_Name\tAsv_Number\tSeq_Number\n")
                if len(lines) > 2:
                    for line in lines[2:]:
                        line = line.strip().split("\t")
                        sample_name = line[0]
                        valied_number = int(line[5])
                        asv_num = sample_number_dict[sample_name]
                        w.write("{}\t{}\t{}\n".format(sample_name, asv_num, valied_number))
        else:
            denoise_path = self.output_dir + "/Deblur_sequence_info.txt"
            with open(mongo_denoise_stat, 'w') as w:
                w.write("Sample_Name\tAsv_Number\tSeq_Number\n")
                for sample in sample_list:
                        if sample in sample_number_dict:
                            asv_num = sample_number_dict[sample]
                        if sample in sequence_number_dict:
                            valied_number = sequence_number_dict[sample]
                        w.write("{}\t{}\t{}\n".format(sample, asv_num, valied_number))
        if os.path.exists(denoise_path):
            os.remove(denoise_path)
        if os.path.exists(mongo_denoise_stat):
            os.rename(mongo_denoise_stat, denoise_path)

    def run(self):
        """
        运行
        :return:
        """
        super(Qiime2DenoiseModule, self).run()
        if self.option("denoise_method") in ["DADA2"]:
            self.run_data2()
        elif self.option("denoise_method") in ["Deblur"]:
            self.run_deblur()

    def end(self):
        """
        结束
        :return:
        """
        super(Qiime2DenoiseModule, self).end()
