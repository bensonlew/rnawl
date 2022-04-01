# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os,re
import shutil
import gevent
from Bio import SeqIO
import numpy as np
import pandas as pd
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.metaasv.common_function import link_dir,link_file


class DeblurDenoiseSplitModule(Module):
    """
    按照样本进行拆分，Deblur 降噪, 基于训练好的数据进行降噪
    """
    def __init__(self, work_id):
        super(DeblurDenoiseSplitModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},##fastq_dir文件夹
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列,自定义模式上传的fasta文件
            {'name': 'database', 'type': 'string'},  # 数据库选择  用于选择参考数据库
            {"name": "min_sample_number", "type": "int", "default": 0},## 过滤长度的值
            {"name": "coverage", "type": "int", "default": 90},  # 根据90%的水平进行筛选， 理论上选择最小的文件
            {"name": "truc_len", "type": "int", "default": 0},##过滤长度的值，默认不过滤，如果不过滤deblur必须要对齐序列
            {"name": "min_size", "type": "int", "default": 1},
            {"name": "min_reads", "type": "int", "default": 10},
        ]
        self.add_option(options)
        self.fastq_dir = os.path.join(self.work_dir, "fastq_dir")
        self.trim_data = self.add_tool("metaasv.select_trimlength")
        self.file_list = []
        self.run_tools = []
        self.merge_tools = []

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("请输入fastq文件夹！")
        if not self.option("database"):
            raise OptionError("请输入数据库类型！")

    def select_coverage(self):
        """
        根据每个样本的全部数据和占全部数据百分比的长度得到每个样本的trim长度
        :return:
        """
        options = {
            "fastq_dir":self.option("fastq_dir"),
            "coverage": self.option("coverage")
            }
        self.trim_data.set_options(options)
        self.trim_data.on("end", self.run_deblur)
        self.trim_data.run()

    def run_deblur(self):
        """
        运行软件deblur，非qiime2流程
        :return:
        """
        self.logger.info("min_sample_number: {}".format(self.option("truc_len")))
        if self.option("truc_len") == 0:
            result_path = os.path.join(self.trim_data.output_dir, "all_result.xls")
            min_number = self.get_min_seq_length(result_path)
        else:
            min_number = self.option("truc_len")
        if os.path.exists(self.fastq_dir):
            shutil.rmtree(self.fastq_dir)
        os.mkdir(self.fastq_dir)
        input_dir = self.option("fastq_dir").prop["path"]
        listdirs = os.listdir(input_dir)
        for file in listdirs:
            self.deblur = self.add_tool("metaasv.deblur_denoise_single")
            file_path = os.path.join(input_dir, file)
            new_file_path = os.path.join(self.fastq_dir, file)
            link_file(file_path, new_file_path)
            options = {
                "input_fastq": new_file_path,
                "database": self.option("database"),
                # "trim_length": self.option("trim_length"),
            }
            if self.option("ref_fasta").is_set:
                self.data_se = os.path.join(self.convert_file.output_dir, "data_se.qza")
                options["input_qza"] = self.data_se
            if self.option("ref_fasta").is_set:
                options["ref_fasta"] = os.path.join(self.taxon_file.output_dir, "ref_fasta.qza")## 用参考的序列
            if self.option("truc_len") != 0:
                options["truc_len"] = self.option("truc_len")
            elif self.option("coverage") != 0:
                options["truc_len"] = min_number
            else:
                options["truc_len"] = 90
            if self.option("min_size"):
                options["min_size"] = self.option("min_size")
            if self.option("min_reads"):
                options["min_reads"] = self.option("min_reads")
            self.deblur.set_options(options)
            self.run_tools.append(self.deblur)
        if len(self.run_tools) > 1:
            self.on_rely(self.run_tools, self.run_merge_file)
        else:
            self.run_tools[0].on('end', self.run_merge_file)
        for tool in self.run_tools:
            tool.run()
            gevent.sleep(0)

    def run_merge_file(self):
        """
        分别对feature-table、fasta和stat统计文件进行合并,fa去重
        :return:
        """
        asv_table = os.path.join(self.work_dir, "asv_table")
        if os.path.exists(asv_table):
            shutil.rmtree(asv_table)
        os.mkdir(asv_table)
        asv_fasta = os.path.join(self.work_dir, "asv_fasta")
        if os.path.exists(asv_fasta):
            shutil.rmtree(asv_fasta)
        os.mkdir(asv_fasta)
        n = 1
        j = 1
        for tool in self.run_tools:
            list_dirs = sorted(os.listdir(tool.output_dir))
            for file in list_dirs:
                file_path = os.path.join(tool.output_dir, file)
                if re.search(r'ASV_reps\.fa', file):
                    new_file = os.path.join(asv_fasta, file + str(n))
                    link_file(file_path, new_file)
                    n += 1
                elif re.search(r'ASV_table\.xls', file):
                    new_file = os.path.join(asv_table, file + str(j))
                    link_file(file_path, new_file)
                    j += 1

        self.merge1 = self.add_tool("metaasv.merge_file")
        options = {
            "input_dir": asv_table,
            "type": "table"
            }
        self.merge1.set_options(options)
        self.merge_tools.append(self.merge1)
        self.merge2 = self.add_tool("metaasv.merge_file")
        options = {
            "input_dir": asv_fasta,
            "type": "fasta"
            }
        self.merge2.set_options(options)
        self.merge_tools.append(self.merge2)

        if len(self.merge_tools) > 1:
            self.on_rely(self.merge_tools, self.create_mongo_data)
        else:
            self.merge_tools[0].on('end', self.create_mongo_data)
        for tool in self.merge_tools:
            tool.run()
            gevent.sleep(1)

    def create_mongo_data(self):
        """
        将ASV_table.xls和降噪后的序列文件统计信息进行合并
        :return:
        """
        result_dada = os.path.join(self.work_dir, "result_dada")
        if os.path.exists(result_dada):
            shutil.rmtree(result_dada)
        os.mkdir(result_dada)
        for tool in self.merge_tools:
            for file in sorted(os.listdir(tool.output_dir)):
                file_path = os.path.join(tool.output_dir, file)
                new_file = os.path.join(result_dada, file)
                link_file(file_path, new_file)
        asv_fasta = os.path.join(result_dada, "ASV_reps.fasta")
        asv_table_origin = os.path.join(result_dada, "ASV_table.xls")
        asv_table = pd.read_table(asv_table_origin, sep="\t", header=0)

        new_gene = {}
        number = 1
        ### 这里先对asv表的丰度进行过滤，然后根据过滤后的asv再提取序列，替换名称
        with open(asv_table_origin, 'r') as f, open(os.path.join(self.output_dir, "ASV_table.xls"), 'w') as w, open(os.path.join(self.output_dir, "ASV_md5.xls"), 'w') as outf:
            outf.write("ASV ID\tmd5\n")
            for line in f:
                if line.startswith("#OTU ID"):
                    line = line.strip().split("\t")
                    line[0] = "ASV ID"
                    new_line = "\t".join(line)
                    w.write(new_line+"\n")
                else:
                    line = line.strip().split("\t")
                    asv_id = line[0]
                    gene_name = "ASV" + str(number)
                    if asv_id not in new_gene:
                        new_gene[asv_id] = gene_name
                    line[0] = gene_name
                    asv_list = [float(x) for x in line[1:]]
                    total_asv_size = int(float(np.sum(asv_list)))
                    # choose_list = filter(lambda x: x != 0.0, asv_list)
                    # if len(choose_list) >= 2: ###在这里筛选asv_size，每条asv在至少两个样本中存在才保留，否则过滤掉
                    # if 10 >= self.option("min_reads"):
                    #     min_reads = 10
                    # else:
                    #     min_reads = self.option("min_reads")
                    if total_asv_size >= int(self.option("min_reads")):## 在这里筛选asv_abundance，所有样本的总和大于等于10的asv才会保留，否则过滤
                        number += 1
                        w.write("\t".join(line) + "\n")
                        outf.write("{}\t{}\n".format(gene_name, asv_id))

        with open(os.path.join(self.output_dir, "ASV_reps.fasta"), "w") as wf:
            for seq_record in SeqIO.parse(asv_fasta, 'fasta'):
                gene_id = seq_record.id
                gene_seq = seq_record.seq
                if gene_id in new_gene:
                    gene_name = new_gene[gene_id]
                    wf.write(">{}\n{}\n".format(gene_name, gene_seq))

        mongo_denoise_stat = self.work_dir + "/denoise_stat.xls"
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

        denoise_path = result_dada + "/Deblur_sequence_info.txt"
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
        self.run_convert_qza()

    def run_convert_qza(self):
        """
        转table转为qza文件
        :return:
        """
        qza_table = os.path.join(self.work_dir, "result_dada")
        allfiles = sorted(os.listdir(qza_table))
        for file in allfiles:
            file_path = os.path.join(qza_table, file)
            self.convert_qza = self.add_tool("metaasv.file_to_qzaqzv")
            if re.search(r'ASV_table', file):
                options = {
                "input_file": file_path,
                "type": "table",
                "prefix": "ASV_table"
                }
            elif re.search(r'ASV_reps', file):
                options = {
                "input_file": file_path,
                "type": "fasta",
                "prefix": "ASV_reps"
                }
            else:
                options = {
                    "input_file": file_path,
                    "type": "table",
                    "prefix": "Deblur_stats"
                }
            self.convert_qza.set_options(options)
            self.file_list.append(self.convert_qza)
        if len(self.file_list) > 1:
            self.on_rely(self.file_list, self.set_output)
        else:
            self.file_list[0].on("end", self.set_output)
        for tool in self.file_list:
            tool.run()
            gevent.sleep(0)

    def get_min_seq_length(self, inputfilr):
        """
        根据传入的inputfile获得最小值
        :param inputfilr:
        :return:
        """
        min_number_list = []
        with open(inputfilr, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                sample_min = int(float(line[1]))
                if sample_min not in min_number_list:
                    min_number_list.append(sample_min)
        min_number = min(min_number_list)
        return min_number

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        dada2_stat = os.path.join(self.output_dir, "Deblur_sequence_info.txt")
        if not os.path.exists(dada2_stat):
            link_file(os.path.join(self.work_dir,"result_dada", "Deblur_sequence_info.txt"),dada2_stat)
        for tool in self.file_list:
            for file in os.listdir(tool.output_dir):
                file_path = os.path.join(tool.output_dir,file)
                out_path = os.path.join(self.output_dir, file)
                link_file(file_path, out_path)
        # for module in self.merge_tools:
        #     for file in os.listdir(module.output_dir):
        #         file_path = os.path.join(module.output_dir,file)
        #         out_path = os.path.join(self.output_dir, file)
        #         link_file(file_path, out_path)
        self.logger.info("设置结果文件目录成功！")
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        super(DeblurDenoiseSplitModule, self).run()
        if self.option("truc_len") == 0:
            self.select_coverage()
        else:
            self.run_deblur()

    def end(self):
        """
        结束
        :return:
        """
        super(DeblurDenoiseSplitModule, self).end()
