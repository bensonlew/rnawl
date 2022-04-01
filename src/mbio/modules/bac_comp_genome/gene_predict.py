#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == "qingchen.zhang"

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import json
import shutil
import random
from mbio.packages.metagbin.common_function import link_dir


class GenePredictModule(Module):
    """
    微生物比较基因组所有样品的基因预测模块
    """
    def __init__(self, work_id):
        super(GenePredictModule, self).__init__(work_id)
        options = [
            {"name": "raw_dir", "type": "string"},###前面已经做过file检查
            {"name": "total_genome", "type": "infile", "format": "sequence.profile_table"},
            {"name": "sequence_dir", "type": "string"}, ##seq序列文件的文件夹
            {"name": "gff_dir", "type": "string"}, ##gff文件的dir文件
            ## json格式：参数比较多按照基因组排序（染色体和质粒逗号分隔），样品有重复：
            # assembly+type+assembly_type+genome_id+genome_type+seq+gff+trna_num+rrna_num+species+gene_prefix
        ]
        self.add_option(options)
        self.run_tools = []

    def check_options(self):
        """
        检查参数
        """
        self.logger.info("正在进行参数检查")
        if not self.option("total_genome"):
            raise OptionError("必须输入流程运行参数文件！")

    def get_params(self):
        """
        根据传入参数分开获得内容
        :return:
        """
        self.get_params_from_list()
        self.majorbio = {}
        self.seq = {}
        self.gff = {}
        self.sample_list = []
        #total_genome = json.loads(self.option("total_genome"))
        #assembly_type = json.loads(self.option("assembly_type"))
        ##确定了参数的顺序格式
        ##assembly+type+assembly_type+genome_id+genome_type+seq+gff+trna_num+rrna_num+species+gene_prefix
        ###共11列
        total_genome_path = self.option("total_genome").prop["path"]
        if self.option("raw_dir"):
            total_path = self.option("raw_dir")
        if self.option("sequence_dir"):
            seq_path = self.option("sequence_dir")
        if self.option("gff_dir"):
            gff_path = self.option("gff_dir")
        with open(total_genome_path, "r") as f:
            for line in f:
                genome = {}
                line = line.strip().split("\t")
                assembly = line[0]
                assembly_type = line[2]
                if assembly_type in ['draft', "DRAFT", "Draft"]: ##如果是扫描图走raw_dir，因为草图不需要拆分序列
                    sample_name = assembly + '.fna'
                    sample_gff_name = assembly + '.gff'
                    if self.option("raw_dir"):
                        new_seq_path = os.path.join(total_path, sample_name)
                        genome["seq"] = new_seq_path
                    else:
                        genome["seq"] = line[5]
                    if self.option("gff_dir"):
                        all_gff = gff_path + assembly
                        new_gff_path = os.path.join(all_gff, sample_gff_name)
                        genome["gff"] = new_gff_path
                    else:
                        genome["gff"] = line[6]
                else:## 如果是完成图走sequence_dir因为完成图是拆分序列的
                    seq_list = []
                    gff_list = []
                    genome_list = line[3].split(",")
                    for ge in genome_list:
                        sample_name = ge + '.fna'
                        sample_gff_name = ge + '.gff'
                        if self.option("sequence_dir"):
                            all_seq = seq_path + assembly
                            new_seq_path = os.path.join(all_seq, sample_name)
                            seq_list.append(new_seq_path)
                        if self.option("gff_dir"):
                            all_gff = gff_path + assembly
                            new_gff_path = os.path.join(all_gff, sample_gff_name)
                            gff_list.append(new_gff_path)
                    if self.option("sequence_dir"):
                        all_seq_path = ','.join(seq_list)
                    else:
                        seq_path = line[5]
                        seq_list.append(seq_path)
                        all_seq_path = ','.join(seq_list)
                    if self.option("gff_dir"):
                        all_gff_path = ','.join(gff_list)
                    else:
                        gff_path = line[6]
                        gff_list.append(gff_path)
                        all_gff_path = ','.join(gff_list)
                    self.logger.info(assembly)
                    genome["seq"] = all_seq_path
                    genome["gff"] = all_gff_path
                genome_id = line[3]
                type = line[1]
                #genome_status = line[3]
                if type in ["seq"]:
                    genome["assembly"] = line[0]
                    genome["assembly_type"] = line[2]
                    genome["genome"] = line[3]
                    genome["genome_type"] = line[4]
                    #genome["seq"] = line[5]
                    genome["species"] = line[-2]
                    genome["gene_tag"] = line[-1]
                    genome["type"] = type
                    self.seq[assembly] = genome
                elif type in ["gff"]:
                    genome["assembly"] = line[0]
                    genome["assembly_type"] = line[2]
                    genome["genome"] = line[3]
                    genome["genome_type"] = line[4]
                    #genome["seq"] = line[5]
                    #genome["gff"] = line[6]
                    genome["trna"] = line[7]
                    genome["rrna"] = line[8]
                    genome["species"] = line[-2]
                    genome["gene_tag"] = line[-1]
                    genome["type"] = type
                    self.gff[assembly] = genome
                elif type in ["majorbio"]:
                    genome["assembly"] = assembly
                    genome["genome"] = genome_id
                    genome["assembly_type"] = assembly_type
                    #genome["genome_status"] = genome_status
                    genome["type"] = type
                    self.majorbio[assembly] = genome
                else:
                    pass
                self.sample_list.append(assembly)
                self.sample_list.sort()


    def get_params_from_list(self):
        """
        从文件夹中的list文件中整理params
        :return:
        """
        pass

    def run_majorbio(self):
        """
        运行从majorbio数据库中获取数据
        :return:
        """
        self.logger.info("正在从数据库中进行单样品的基因预测")
        for sample in self.sample_list:
            if sample in self.majorbio.keys():
                genome_dict = {}
                params = self.majorbio[sample]
                predict_seq = self.add_module("bac_comp_genome.gene_predict_single") #单个样品的基因预测，非编码预测，重复序列预测
                assembly = params["assembly"]
                assembly_type = params["assembly_type"]
                genomes = params["genome"]
                genome_dict[assembly] = genomes
                sample = json.dumps(genome_dict)
                opts = {
                    "sample": sample,
                    "type": "majorbio",
                    "assembly_type": assembly_type,
                    }
                predict_seq.set_options(opts)
                self.run_tools.append(predict_seq)

    def run_seq(self):
        """
        根据序列seq直接分析获取数据
        :return:
        """
        self.logger.info("正在对上传序列进行基因预测")
        for sample in self.sample_list:
            if sample in self.seq.keys():
                params = self.seq[sample]
                genome_dict = {}
                predict_seq = self.add_module("bac_comp_genome.gene_predict_single")  # 单个样品的基因预测，非编码预测，重复序列预测
                assembly = params["assembly"]
                assembly_type = params["assembly_type"]
                genomes = params["genome"]
                genomes_list = params["genome"].split(",")
                genome_type = params["genome_type"].split(",")
                seqs = params["seq"].split(",")
                seq_dict = dict(zip(genomes_list, seqs))
                genome_dict[assembly] = genomes
                if params["gene_tag"] != "-":
                    gene_tag = params["gene_tag"].split(",")
                else:
                    gene_tag = params["gene_tag"]
                gene_dict = dict(zip(genomes_list, gene_tag))
                gene_type_dict = dict(zip(genomes_list, genome_type))
                type = params["type"]
                self.logger.info("aaaa++++++++++++++++++++%s"%gene_dict)
                opts = {
                    "sample": json.dumps(genome_dict),
                    "genome_type": json.dumps(gene_type_dict),
                    "assembly_type": assembly_type,
                    "genome_prefix": json.dumps(gene_dict),
                    "seq": json.dumps(seq_dict),
                    "type": type
                }
                predict_seq.set_options(opts)
                self.run_tools.append(predict_seq)

    def run_gff(self):
        """
        根据gff+seq序列获取数据
        :return:
        """
        self.logger.info("正在对上传的gff文件进行基因预测")
        for sample in self.sample_list:
            if sample in self.gff.keys():
                params = self.gff[sample]
                genome_dict = {}
                predict_seq = self.add_module("bac_comp_genome.gene_predict_single") #单个样品的基因预测，非编码预测，重复序列预测
                assembly = params["assembly"]
                assembly_type = params["assembly_type"]
                genomes = params["genome"]
                genomes_list = params["genome"].split(",")
                genome_type = params["genome_type"].split(",")
                gff_all =  params["gff"].split(",")
                gff_dict = dict(zip(genomes_list, gff_all))
                seqs = params["seq"].split(",")
                seq_dict = dict(zip(genomes_list, seqs))
                genome_dict[assembly] = genomes
                if params["gene_tag"] != "-":
                    gene_tag = params["gene_tag"].split(",")
                else:
                    gene_tag = params["gene_tag"]
                gene_dict = dict(zip(genomes_list, gene_tag))
                gene_type_dict = dict(zip(genomes_list, genome_type))
                trna = params["trna"]
                rrna = params["rrna"]
                type = params["type"]
                opts = {
                    "sample": json.dumps(genome_dict),
                    "genome_type": json.dumps(gene_type_dict),
                    "assembly_type": assembly_type,
                    "genome_prefix": json.dumps(gene_dict),
                    #"fasta_dir": self.option("raw_dir"),
                    "seq": json.dumps(seq_dict),
                    "gff": json.dumps(gff_dict),
                    "trna_num": trna,
                    "rrna_num": rrna,
                    'type': type
                }
                predict_seq.set_options(opts)
                self.run_tools.append(predict_seq)

    def run_seq_stat(self):
        """
        对seq和gff两种情况的基因组进行统计基础信息
        :return:
        """
        if len(self.gff) != 0 or len(self.seq) != 0:
            self.seq_stat = self.add_tool("bac_comp_genome.seq_stat")
            self.logger.info("正在对序列进行统计")
            opts = {
                "total_genome": self.option("total_genome"),
                }
            if self.option("raw_dir"):
                total_path = self.option("raw_dir")
                opts["raw_dir"] = total_path
            if self.option("sequence_dir"):
                seq_path = self.option("sequence_dir")
                opts["sequence_dir"] = seq_path
            self.seq_stat.set_options(opts)
            self.seq_stat.on("end", self.set_output)
            self.seq_stat.run()
        else:
            self.set_output()

    def run(self):
        self.logger.info("正在开始进行基因预测")
        super(GenePredictModule, self).run()
        self.get_params()
        #if len(self.sample_list) == 0:
           # self.set_error("所传的样本list为空")
        if len(self.majorbio) != 0:
            self.run_majorbio()
            self.logger.info("开始从数据库中运行啦")
        if len(self.gff) != 0:
            self.run_gff()
        if len(self.seq) != 0:
            self.run_seq()
        self.on_rely(self.run_tools, self.run_seq_stat)
        for module in self.run_tools:
            module.run()
            gevent.sleep(0)


    def write_list(self):
        """
        output结果文件下生成文件夹的list，以确保后面的样本文件夹的顺序是一样的
        :return:
        """
        file_path = self.output_dir + "/list.txt"
        with open(file_path, "w") as f:
            f.write("Sample_name\tDir_path\n")
            for sample_name in self.sample_list:
                dir_path = os.path.join(self.output_dir, sample_name)
                f.write("{}\t{}\n".format(sample_name, dir_path))

    def set_output(self):
        """
        将各个模块的结果输出至output
        """
        for i in self.run_tools:
            for j in os.listdir(i.output_dir):
                old_dir = os.path.join(i.output_dir, j)
                new_dir = os.path.join(self.output_dir, j)
                if os.path.exists(new_dir):
                    shutil.rmtree(new_dir)
                link_dir(old_dir, new_dir)
        if os.path.exists(self.output_dir + "list.txt"):
            os.remove(self.output_dir + "list.txt")
        # self.write_list()
        if os.path.exists(os.path.join(self.output_dir, "upload_file.xls")):
            os.remove(os.path.join(self.output_dir, "upload_file.xls"))
        if len(self.gff) != 0 or len(self.seq) != 0:
            os.link(os.path.join(self.seq_stat.output_dir, "upload_file.xls"), os.path.join(self.output_dir, "upload_file.xls"))
        self.end()

    def end(self):
        super(GenePredictModule, self).end()

    def get_list(self):
        """
        从样品list文件中获得样本顺序
        :return:
        """
        sample_dict = json.loads(self.option("assembly_type"))
        self.sample_list = sorted(sample_dict)

    def get_gff(self):
        """
        从remote_dir文件夹下将list.txt获取样品
        :return:
        """
        sp_list = []
        seq_dict = {}
        gff_dict = {}
        dir_path = self.option("raw_dir").prop["path"]
        list_path = os.path.join(dir_path, "list.txt")
        with open(list_path, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                sample = line[0]
                file = line[1]
                if re.search(r"\.fa", line[1]):
                    seq_path = os.path.join(dir_path, line[1])
                elif re.search(r"\.gff", line[1]):
                    gff_path = os.path.join(dir_path, line[1])
                type = line[2]
                genome_id = line[3]
                if sample not in sp_list:
                    sp_list.append(sample)
                if genome_id not in seq_dict.keys():
                    seq_dict[genome_id] = seq_path
                if genome_id not in gff_path.keys():
                    gff_dict[genome_id] = gff_path


    def get_sequence_type(self,file):
        dict ={}
        type =[]
        with open(file, "rb") as l:
            raw_lines = l.readlines()
            for line in raw_lines[1:]:
                line2 = line.strip("\r\n").split("\t")
                dict[line2[5]] = line2[5]
        for k in dict.iterkeys():
            type.append(k)
        if len(type) ==1:
            sequence_type=type[0]
        else:
            sequence_type = ",".join(type)
        return sequence_type

    def get_gene_prefix(self):
        de =""
        prefix ={}
        if self.option("analysis") in ["complete"]:
            file =self.assemble_assess.work_dir + "/plasmid.type.xls"
            with open(file,"r") as f:
                lines =f.readlines()
                for line in lines[0:]:
                    line =line.rstrip("\r\n").split("\t")
                    prefix[line[0]]=line[1]
            f.close()
            de =str(prefix)
        return de

    def get_seq_type(self,file_dir):
        list =[]
        files =os.listdir(file_dir)
        for file in files:
            lst=file.split(".")
            list.append(lst[0])
        if len(list) == 1:
            seq_type =list[0]
        else:
            seq_type = ",".join(list)
        return seq_type