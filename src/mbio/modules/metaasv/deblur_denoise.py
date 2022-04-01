# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os,re
import shutil
import gevent
from Bio import SeqIO
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.metaasv.common_function import link_dir,link_file


class DeblurDenoiseModule(Module):
    """
    Deblur 降噪, 基于训练好的数据进行降噪
    """
    def __init__(self, work_id):
        super(DeblurDenoiseModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},##fastq_dir文件夹
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列,自定义模式上传的fasta文件
            {'name': 'database', 'type': 'string'},  # 数据库选择  用于选择参考数据库
            {"name": "min_sample_number", "type": "int", "default": 0},## 样本的最小样本序列数
        ]
        self.add_option(options)
        # self.config = os.path.join(Config().SOFTWARE_DIR, "database/taxon_db/qiime2_qza")
        self.file_list = []

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("请输入fastq文件夹！")
        if not self.option("database"):
            raise OptionError("请输入数据库类型！")

    def run_convert_file(self):
        """
        将fastq_dir文件夹转为qza文件
        :return:
        """
        self.convert_file = self.add_tool("metaasv.file_to_qza")
        options = {
            "input_file":self.option("fastq_dir").prop['path'],
            "type": "fq_dir",
            "prefix": "data_se"
            }
        self.convert_file.set_options(options)
        if self.option("ref_fasta").is_set:
            self.convert_file.on("end", self.run_taxon_file)
        else:
            self.convert_file.on("end", self.run_deblur)
        self.convert_file.run()

    def run_taxon_file(self):
        """
        将ref_fasta文件转为qza文件
        :return:
        """
        self.taxon_file = self.add_tool("metaasv.file_to_qza")
        options = {
            "input_file":self.option("ref_fasta").prop['path'],
            "type": "fasta",
            "prefix": "ref_fasta"
            }
        self.taxon_file.set_options(options)
        self.taxon_file.on("end", self.run_deblur)
        self.taxon_file.run()

    def run_deblur(self):
        """
        运行软件qiime2 deblur_denoise
        :return:
        """
        self.deblur = self.add_tool("metaasv.deblur_denoise")
        options = {
            "fastq_dir": self.option("fastq_dir"),
            "database": self.option("database")
            }
        if self.option("ref_fasta").is_set:
            self.data_se = os.path.join(self.convert_file.output_dir, "data_se.qza")
            options["input_qza"] = self.data_se
        if self.option("ref_fasta").is_set:
            options["ref_fasta"] = os.path.join(self.taxon_file.output_dir, "ref_fasta.qza")## 用参考的序列
        # else:
        #     ref_fasta = os.path.join(self.config, self.option("database"), "ref-seqs.qza")## 用原始数据库的序列
        #     options["ref_fasta"] = ref_fasta
        if self.option("min_sample_number") != 0:
            options["min_sample_number"] = self.option("min_sample_number")
        self.deblur.set_options(options)
        self.deblur.on("end", self.run_convert_qza)
        self.deblur.run()

    def run_convert_qza(self):
        """
        转qza文件转为table
        :return:
        """
        qza_table = os.path.join(self.work_dir, "qza_table")
        if os.path.exists(qza_table):
            shutil.rmtree(qza_table)
        os.mkdir(qza_table)
        for file in os.listdir(self.deblur.output_dir):
            if re.search(r".qza", file):
                file_path = os.path.join(self.deblur.output_dir, file)
                new_file_path = os.path.join(qza_table, file)
                link_file(file_path, new_file_path)
        allfiles = os.listdir(qza_table)
        for file in allfiles:
            file_path = os.path.join(qza_table, file)
            self.convert_qza = self.add_tool("metaasv.qza_to_table")
            if re.search(r'ASV_table', file):
                options = {
                "input_file": file_path,
                "type": "qza_biom",
                "prefix": file.strip(".qza")
                }
            else:
                options = {
                    "input_file": file_path,
                    "type": "qza",
                    "prefix": file.strip(".qza")
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

    def replace_name(self, asv_fasta, asv_table):
        """
        根据asv_fasta和ASV_table，生成
        :param asv_fasta:
        :param asv_table:
        :return:
        """
        origin_gene_dict = {}
        new_gene = {}
        number = 1
        with open(os.path.join(self.output_dir, "ASV_reps2.fasta"), "w") as wf, open(os.path.join(self.output_dir, "ASV_md5.xls"), "w") as outf:
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

        with open(asv_table, 'r') as f, open(os.path.join(self.output_dir, "ASV_table2.xls"), 'w') as w:
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
                        w.write("\t".join(line) + "\n")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        link_dir(self.deblur.output_dir, self.output_dir)
        for tool in self.file_list:
            for file in os.listdir(tool.output_dir):
                file_path = os.path.join(tool.output_dir,file)
                out_path = os.path.join(self.output_dir, file)
                link_file(file_path, out_path)
        self.logger.info("设置结果文件目录！")
        asv_fasta = os.path.join(self.output_dir, "ASV_reps.fasta")
        asv_table = os.path.join(self.output_dir, "ASV_table.xls")
        if os.path.exists(os.path.join(self.output_dir, "dna-sequences.fasta")):
            os.rename(os.path.join(self.output_dir, "dna-sequences.fasta"), os.path.join(self.output_dir, "ASV_reps.fasta"))
        if os.path.exists(os.path.join(self.output_dir, "feature-table.biom")):
            os.rename(os.path.join(self.output_dir, "feature-table.biom"), os.path.join(self.output_dir, "ASV_table.biom"))
        if os.path.exists(os.path.join(self.output_dir, "stats.tsv")):
            os.rename(os.path.join(self.output_dir, "stats.tsv"), os.path.join(self.output_dir, "Deblur_sequence_info.txt"))
        self.replace_name(asv_fasta, asv_table)
        if os.path.exists(os.path.join(self.output_dir, "ASV_reps2.fasta")):
            os.rename(os.path.join(self.output_dir, "ASV_reps2.fasta"), asv_fasta)
        if os.path.exists(os.path.join(self.output_dir, "ASV_table2.xls")):
            os.rename(os.path.join(self.output_dir, "ASV_table2.xls"), asv_table)
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        super(DeblurDenoiseModule, self).run()
        if self.option("ref_fasta").is_set:
            self.run_taxon_file()
        else:
            self.run_deblur()

    def end(self):
        """
        结束
        :return:
        """
        super(DeblurDenoiseModule, self).end()
