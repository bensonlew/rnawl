# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os,re
import shutil
import gevent
import numpy as np
from Bio import SeqIO
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file


class Dada2DenoiseModule(Module):
    """
    按照样本并行 运行dada2软件 降噪
    """
    def __init__(self, work_id):
        super(Dada2DenoiseModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},##fastq_dir文件夹
            {"name": "max_ee", "type": "int", "default": 2},
            {"name": "trunc_q", "type": "int", "default": 0}, ##判断打断序列的大小，如果为0则不打断，即不质控
        ]
        self.add_option(options)
        self.fastq_dir = os.path.join(self.work_dir, "fastq_dir")
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

    def run_dada2(self):
        """
        运行软件qiime2 data2_denoise
        :return:
        """
        if os.path.exists(self.fastq_dir):
            shutil.rmtree(self.fastq_dir)
        os.mkdir(self.fastq_dir)
        input_dir = self.option("fastq_dir").prop["path"]
        listdirs = os.listdir(input_dir)
        for file in listdirs:
            self.dada2 = self.add_tool("metaasv.dada2_denoise")
            file_path = os.path.join(input_dir, file)
            new_file_path = os.path.join(self.fastq_dir, file)
            link_file(file_path, new_file_path)
            options = {
                "input_fastq": new_file_path,
                "cpu": 3,
                "max_ee": self.option("max_ee"),
                "trunc_q": self.option("trunc_q"),
                }
            self.dada2.set_options(options)
            self.run_tools.append(self.dada2)
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
        asv_stat = os.path.join(self.work_dir, "asv_stat")
        if os.path.exists(asv_stat):
            shutil.rmtree(asv_stat)
        os.mkdir(asv_stat)
        n = 1
        j = 1
        i = 1
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
                elif re.search(r'DADA2_stats\.xls', file):
                    new_file = os.path.join(asv_stat, file + str(i))
                    link_file(file_path, new_file)
                    i += 1
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
        self.merge3 = self.add_tool("metaasv.merge_file")
        options = {
            "input_dir": asv_stat,
            "type": "stat"
            }
        self.merge3.set_options(options)
        self.merge_tools.append(self.merge3)
        if len(self.merge_tools) > 1:
            self.on_rely(self.merge_tools, self.replace_name)
        else:
            self.merge_tools[0].on('end', self.replace_name)
        for tool in self.merge_tools:
            tool.run()
            gevent.sleep(1)

    def replace_name(self):
        """
        根据asv_fasta和ASV_table，生成最终的结果表,并生成md5的结果表
        :return:
        :param filter_size: 过滤asv_size的大小，低于这个值的asv将会被过滤掉
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
        asv_table = os.path.join(result_dada, "ASV_table.xls")
        origin_gene_dict = {}
        new_gene = {}
        number = 1
        with open(os.path.join(self.output_dir, "ASV_reps.fasta"), "w") as wf, open(os.path.join(self.output_dir, "ASV_md5.xls"), 'w') as outf:
            outf.write("ASV ID\tmd5\n")
            for seq_record in SeqIO.parse(asv_fasta, 'fasta'):
                gene_id = seq_record.id
                gene_seq = seq_record.seq
                origin_gene_dict[gene_id] = gene_seq
                gene_name = "ASV" + str(number)
                new_gene[gene_id] = gene_name
                wf.write(">{}\n{}\n".format(gene_name, gene_seq))
                number += 1
                outf.write("{}\t{}\n".format(gene_name, gene_id))

        with open(asv_table, 'r') as f, open(os.path.join(self.output_dir, "ASV_table.xls"), 'w') as w:
            for line in f:
                if line.startswith("#OTU ID"):
                    line = line.strip().split("\t")
                    line[0] = "ASV ID"
                    new_line = "\t".join(line)
                    w.write(new_line+"\n")
                else:
                    line = line.strip().split("\t")
                    asv_id = line[0]
                    if asv_id in new_gene.keys():
                        gene_name = new_gene[asv_id]
                        line[0] = gene_name
                        asv_list = [float(x) for x in line[1:]]
                        total_asv_size = int(float(np.sum(asv_list)))
                        if total_asv_size >= 2:
                            w.write("\t".join(line) + "\n")
        self.run_convert_qza()

    def run_convert_qza(self):
        """
        转table转为qza文件
        :return:
        """
        qza_table = os.path.join(self.work_dir, "qza_table")
        if os.path.exists(qza_table):
            shutil.rmtree(qza_table)
        os.mkdir(qza_table)
        link_file(os.path.join(self.work_dir, "result_dada", "ASV_table.xls"), os.path.join(self.work_dir, "qza_table", "ASV_table.xls"))
        link_file(os.path.join(self.work_dir, "result_dada", "ASV_reps.fasta"), os.path.join(self.work_dir, "qza_table", "ASV_reps.fasta"))
        link_file(os.path.join(self.work_dir, "result_dada", "DADA2_sequence_info.txt"), os.path.join(self.work_dir, "qza_table", "DADA2_sequence_info.txt"))
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
                    "prefix": "DADA2_stats"
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

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        dada2_stat = os.path.join(self.output_dir, "DADA2_sequence_info.txt")
        if not os.path.exists(dada2_stat):
            link_file(os.path.join(self.work_dir,"qza_table", "DADA2_sequence_info.txt"),dada2_stat)
        # link_dir(self.data2.output_dir, self.output_dir)
        for tool in self.file_list:
            for file in os.listdir(tool.output_dir):
                file_path = os.path.join(tool.output_dir,file)
                out_path = os.path.join(self.output_dir, file)
                link_file(file_path, out_path)
        self.logger.info("设置结果文件目录成功！")
        # asv_fasta = os.path.join(self.output_dir, "ASV_reps.fasta")
        # asv_table = os.path.join(self.output_dir, "ASV_table.xls")
        # os.rename(os.path.join(self.output_dir, "dna-sequences.fasta"), os.path.join(self.output_dir, "ASV_reps.fasta"))
        # os.rename(os.path.join(self.output_dir, "feature-table.biom"), os.path.join(self.output_dir, "ASV_table.biom"))
        # os.rename(os.path.join(self.output_dir, "stats.tsv"), os.path.join(self.output_dir, "DADA2_sequence_info.txt"))
        # self.replace_name(asv_fasta, asv_table)
        # os.rename(os.path.join(self.output_dir, "ASV_reps2.fasta"), asv_fasta)
        # os.rename(os.path.join(self.output_dir, "ASV_table2.xls"), asv_table)
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        super(Dada2DenoiseModule, self).run()
        self.run_dada2()

    def end(self):
        """
        结束
        :return:
        """
        super(Dada2DenoiseModule, self).end()
