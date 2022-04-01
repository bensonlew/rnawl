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
# from mbio.packages.bac_comp_genome.get_mongo_stat import get_stat_from_mongo



class GenePredictSingleModule(Module):
    """
    微生物比较基因组单样品基因预测模块
    1. 实现单样品3种情况的预测：majorbio、seq、gff；
    2. 3种情况并行运行，并当为完成图的情况时染色体和质粒是分开并行跑的；
    3. 并行跑完进行合并，合并完成进行基因组水平上的统计
    """
    def __init__(self, work_id):
        super(GenePredictSingleModule, self).__init__(work_id)
        options = [
            {"name": "sample", "type": "string"},###样本名称及对应的genome{"GCF_000217215.1": "NC_005823.1", "GCF_000217215.1": "NC_005867.1"}
            {"name": "genome_prefix", "type": "string"},  ###基因组的前缀，比如{"NC_005823.1":"LEP2GSC107_RS", "NC_005867.1":"LEP2GSC107_RS"}
            {"name": "genome_type", "type": "string"},###完成图基因组的类型，是chromosome还是plasmid，扫描图可以不同此字段；形式{"NC_005823.1": "chromosome", "NC_005867.1": "plasmid"}
            {"name": "assembly_type", "type": "string", "default": "draft"},###组装序列是完成图还是扫描图{}
            {"name": "trna_num", "type": "string", "default": ""}, ##上传的gff文件中是否有trna；
            {"name": "rrna_num", "type": "string", "default": ""}, ##上传的gff文件中是否有rRNA；
            {"name": "seq", "type": "string"}, ###{"NC_005823.1":"path1", "NC_005867.1":"path2"}
            {"name": "gff", "type": "string"}, ###{"NC_005823.1":"path1", "NC_005867.1":"path2"}
            {"name": "type", "type": "string"}, ###是三种情况的哪种情况，majorbio\seq\gff
            #{"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},###fasta文件dir
            {"name": "species_name", "type": "string"}, ###物种名称，如果传过来的基因前缀名称为空，则取物种名称前三位+三位数字
            {"name": "single_genome", "type": "infile", "format":"sequence.profile_table"}, #参数的infile
            {"name": "software_list", "type":"string", "default":"glimmer"},  #使用的软件，逗号分割
            {"name": "p_software_list", "type":"string", "default":"genemark"},  #使用的软件，逗号分割
        ]
        self.add_option(options)
        self.step.add_steps("gene_predict")
        self.tools = []
        self.merge_tools = []

    def check_options(self):
        self.logger.info("正在进行参数检查")
        if not self.option("sample"):
            raise OptionError("必须提供样品名称")
        if not self.option("type"):
            raise OptionError("必须提供seq|gff|majorbio三种方式任何一种")
        else:
            sample = self.get_params()
            if self.option("type") in ["seq"]:
                self.seq = json.loads(self.option("seq"))
                self.gene_tag = json.loads(self.option("genome_prefix"))
                self.genome_type = json.loads(self.option("genome_type"))
                for i in self.genome_list:
                    self.logger.info("")
                    if not os.path.exists(self.seq[i]):
                        raise OptionError("必须提供样品基因组序列")
                    if not self.genome_type.has_key(i):
                        raise OptionError("必须提供样品基因组类型-chromosome")
            elif self.option("type") in ["gff"]:
                self.gene_tag = json.loads(self.option("genome_prefix"))
                self.genome_type = json.loads(self.option("genome_type"))
                self.gff = json.loads(self.option("gff"))
                self.seq = json.loads(self.option("seq"))
                for i in self.genome_list:
                    if not os.path.exists(self.seq[i]):
                        raise OptionError("必须提供样品基因组序列")
            elif self.option("type") in ["majorbio"]:
                #self.genome_status = json.loads(self.option("genome_type"))

                for key in self.genome_list:
                    if not key:
                        raise OptionError("必须提供样品基因组名称")
                    #if not self.genome_type[key]:
                        #raise OptionError("必须提供样品%s基因组%s的类型"%(sample, key))###这里是chromosome还是plasmid

        if not self.option("assembly_type"):###是完成图还是扫描图
            raise OptionError("必须提供样品类型")

    def get_params(self):
        """
        获取参数
        :return:
        """
        if self.option("single_genome").is_set:
            self.logger.info("哈哈哈哈")
        else:
            sample = json.loads(self.option("sample").encode("raw_unicode_escape"))
            self.convert_type()
            sample_name = []
            genome_id =[]
            for key in sample.keys():
                genome_id = sample[key].split(",")
                sample_name.append(key)
            self.genome_list = genome_id
            name = "".join(set(sample_name))
            return name

    def run_majorbio(self):
        """
        并行运行基因组预测
        :return:
        """
        self.convert_type()
        self.logger.info("正在从mongo中获取统计结果")
        #self.get_majorbio_stat()
        self.logger.info("正在从数据库中进行基因预测")
        for genome in self.genome_list:
            self.majorbio = self.add_tool("bac_comp_genome.mongo_gene_predict")
            opts = {
                "sample": self.sample,
                "genome": genome,
                "genome_type": self.assembly_type,
            }
            self.majorbio.set_options(opts)
            self.tools.append(self.majorbio)

    def run_seq(self):
        """
        根据上传序列预测
        :return:
        """
        self.logger.info("正在根据上传序列seq进行基因预测")
        self.seq = json.loads(self.option("seq").encode("raw_unicode_escape"))
        seq_number = 1
        for genome in self.genome_list:
            opts = {}
            self.predict_seq = self.add_module("bac_comp_genome.gene_predict_seq")
            if genome in self.gene_tag.keys():
                gene_tag = self.gene_tag[genome].encode("raw_unicode_escape")
            else:
                seq_number += 1
                if self.option("species_name") != "-":
                    gene_tag = str(self.option("species_name"))[0:2] + str("%03d"%(seq_number))
                else:
                    gene_tag = str(genome)[0:2] + str("%03d"%(seq_number))
            seq = self.seq[genome].encode("raw_unicode_escape")
            self.logger.info("%s" %(seq))
            self.logger.info("%s" %(gene_tag))
            self.logger.info("aaaa+++++++++++%s"%gene_tag)
            if self.assembly_type in ["draft"]:
                opts = {
                    "genome": seq,
                    "genome_name": genome, #genome名称或者id
                    "gene_prefix":gene_tag,
                    "genome_type": self.assembly_type,
                    "sample": self.sample,
                    "software_list": self.option("software_list"),
                }
            elif self.assembly_type in ["complete"]:
                genome_type = self.genome_type[genome]
                if genome_type in["chromosome", "Chromosome", "CHROMOSOME"]:
                    opts = {
                        "genome": seq,
                        "genome_name": genome,#genome名称
                        "gene_prefix":gene_tag,
                        "genome_type": self.assembly_type,
                        "sample": self.sample,
                        "software_list": self.option("software_list"),
                    }
                elif genome_type in["plasmid", "Plasmid", "PLASMID"]:
                    opts = {
                        "genome_plas": seq,
                        "genome_name": genome,#genome名称
                        "plasmid_prefix" : gene_tag,
                        "genome_type": self.assembly_type,
                        "sample": self.sample,
                        "p_software_list": self.option("p_software_list"),
                    }
            self.logger.info(opts)
            self.predict_seq.set_options(opts)
            self.tools.append(self.predict_seq)

    def run_gff(self):
        """
        根据上传的gff+seq进行预测
        :return:
        """
        self.logger.info("正在根据上传的seq+gff文件进行基因预测")
        self.seq = json.loads(self.option("seq"))
        self.gff = json.loads(self.option("gff"))
        seq_number = 1
        for genome in self.genome_list:
            self.predict_gff = self.add_module("bac_comp_genome.gene_predict_gff")
            gff = self.gff[genome]
            seq = self.seq[genome]
            if genome and os.path.exists(gff):
                if genome in self.gene_tag.keys():
                    gene_tag = self.gene_tag[genome]
                else:
                    seq_number += 1
                    if self.option("species_name") != "-":
                        gene_tag = str(self.option("species_name"))[0:2] + str("%03d" % (seq_number))
                    else:
                        gene_tag = str(genome)[0:2] + str("%03d" % (seq_number))
                genome_type = self.assembly_type  # 这里做一个转变只输入基因组的类型就行了
                if self.assembly_type in ["draft"]:
                    opts = {
                        "genome": seq,
                        "input_gff": gff,
                        "sample": self.sample,
                        "gene_prefix": gene_tag,
                        "genome_type": genome_type,
                        "genome_name": genome,
                        "trna": str(self.option("trna_num")),
                        "rrna": str(self.option("rrna_num"))
                    }
                    self.logger.info(opts)
                    self.predict_gff.set_options(opts)
                    self.tools.append(self.predict_gff)
                elif self.assembly_type in ["complete", "chromosome"]:
                    genome_type = self.genome_type[genome]
                    self.logger.info(genome_type)
                    if genome_type in["chromosome", "Chromosome", "CHROMOSOME"]:
                        opts = {
                            "genome": seq,
                            "input_gff": gff,
                            "sample": self.sample,
                            "gene_prefix": gene_tag,
                            "genome_type": genome_type,
                            "genome_name": genome,
                            "trna": str(self.option("trna_num")),
                            "rrna": str(self.option("rrna_num"))
                        }
                    elif genome_type in["plasmid", "Plasmid", "PLASMID"]:
                        opts = {
                            "genome_plas": seq,
                            "input_gff": gff,
                            "genome_name": genome,#genome名称
                            "plasmid_prefix" : gene_tag,
                            "genome_type": self.assembly_type,
                            "sample": self.sample,
                            "trna": str(self.option("trna_num")),
                            "rrna": str(self.option("rrna_num"))
                        }
                    self.logger.info(opts)
                    self.predict_gff.set_options(opts)
                    self.tools.append(self.predict_gff)



    def merge_file(self):
        """
        将完成图的质粒和染色体所有序列合并起来
        :return:
        """
        if os.path.exists(self.work_dir + "/predict_dir"):
            shutil.rmtree(self.work_dir + "/predict_dir")
        os.mkdir(self.work_dir + "/predict_dir")
        for tool in self.tools:
            for file in os.listdir(tool.output_dir):
                old_path = os.path.join(tool.output_dir, file)
                new_path = os.path.join(self.work_dir + "/predict_dir", file)
                if os.path.exists(new_path):
                    os.remove(new_path)
                if re.search(r"rRNA.faa|tRNA.faa", old_path):
                    pass
                else:
                    os.link(old_path, new_path)
        merge_dir = self.work_dir + "/predict_dir"
        gff_list = ["_CDS.gff", "_rRNA.gff","_tRNA.gff"]
        if self.option("type") in ["seq"]:
            gff_list.append("_stat.xls")
        seq_list = ["_CDS.faa", "_CDS.fna", "_16S.fna", "_rRNA.fna", "_tRNA.fna"]
        for file in gff_list:
            self.logger.info(file)
            merge_tool = self.add_tool("bac_comp_genome.cat_table")
            opts = {
                "merge_dir": merge_dir,
                "prefix": self.sample,
                "header": "True",
                "table_name": file,

                }
            merge_tool.set_options(opts)
            self.merge_tools.append(merge_tool)
        for fnn in seq_list:
            merge_tool = self.add_tool("bac_comp_genome.cat_table")
            self.logger.info(fnn)
            opts = {
                "merge_dir": merge_dir,
                "prefix": self.sample,
                "header": "False",
                "table_name": fnn
                }
            merge_tool.set_options(opts)
            self.merge_tools.append(merge_tool)
        if len(self.merge_tools) >= 1:
            self.on_rely(self.merge_tools, self.predict_stat)
        for tool in self.merge_tools:
            tool.run()
            gevent.sleep(0)

    def predict_stat(self):
        """
        对所有的gff文件进行个数统计
        :return:
        """
        merge_dir = ""
        self.logger.info("开始进行majorbiodb基因预测统计")
        if len(self.genome_list) > 1:
            merge_dir = self.get_all_file()
        elif len(self.genome_list) == 1:
            self.logger.info("正在转换文件夹")
            if self.option("type") in ["majorbio"]:
                merge_dir = self.majorbio.output_dir
            elif self.option("type") in ["seq"]:
                merge_dir = self.predict_seq.output_dir
            elif self.option("type") in ["gff"]:
                merge_dir = self.predict_gff.output_dir
            else:
                self.logger.info("运行结果为空")
        else:
            self.set_error("未能获取基因预测的gff文件夹")
        self.stat = self.add_tool("bac_comp_genome.predict_stat")
        opts = {
            "merge_dir": merge_dir,
            "sample": self.sample
            }
        self.stat.set_options(opts)
        self.stat.on("end", self.set_output)
        self.stat.run()

    def convert_type(self):
        """
        需要对组装类型进行转换
        :return:
        """
        if self.option("assembly_type") in ["Draft", "draft", "DRAFT"]:
            self.assembly_type = "draft"
        elif self.option("assembly_type") in ["Complete Genome", "complete", "Complete", "COMPLETE"]:
            self.assembly_type = "complete"
        elif self.option("assembly_type") in ["Chromosome", "chromosome"]:
            self.assembly_type = "chromosome"
        else:
            self.assembly_type = ""

    def get_all_file(self):
        """
        获得前面所有merge后的结果文件
        :return:
        """
        merge_dir = self.work_dir + "/merge_dir"
        if os.path.exists(merge_dir):
            shutil.rmtree(merge_dir)
        os.mkdir(merge_dir)
        for tool in self.merge_tools:
            for file in os.listdir(tool.output_dir):
                old_path = os.path.join(tool.output_dir, file)
                new_file = os.path.join(merge_dir, file)
                if os.path.exists(new_file):
                    os.remove(new_file)
                os.link(old_path, new_file)
        return merge_dir

    def get_sample(self):
        """
        获得样本list
        :return:
        """
        self.sample_list = []
        if self.option("fasta_dir").is_set:
            dir_path = self.option("fasta_dir").prop["path"]
            list_path = os.path.join(dir_path, "/list.txt")
            with open(list_path, "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    self.sample_list.append(line[0])
        else:
            sample_dict = json.loads(self.option("sample"))
            self.sample_list = sorted(sample_dict)


    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在设置结果文件目录")
        sample = self.sample
        if os.path.exists(self.output_dir + "/" + sample):
            shutil.rmtree(self.output_dir + "/" + sample)
        os.mkdir(self.output_dir + "/" + sample)
        result_path = self.output_dir + "/" + sample + "__bak"
        last_path = self.output_dir + "/" + sample
        if len(self.tools) == 1:
            if self.option("type") in ["majorbio"]:
                link_dir(self.majorbio.output_dir, result_path)
            elif self.option("type") in ["seq"]:
                link_dir(self.predict_seq.output_dir, result_path)
            elif self.option("type") in ["gff"]:
                link_dir(self.predict_gff.output_dir, result_path)
            if os.path.exists(last_path):
                shutil.rmtree(last_path)
            os.mkdir(last_path)
            for file in os.listdir(result_path):
                self.logger.info("开始对文件名称进行改名")
                file_path = os.path.join(result_path, file)
                if re.search(r"16S\.fna", file):
                    self.logger.info("%s"%file)
                    file_name = self.sample + "_16S.fna"
                elif re.search(r"sample_stat\.xls", file):
                    file_name = self.sample + "_sample_stat.xls"
                elif re.search(r"_stat\.xls", file):
                    file_name = self.sample + "_sample_stat.xls"
                elif re.search(r"CDS\.faa", file):
                    file_name = self.sample + "_CDS.faa"
                elif re.search(r"CDS\.fna", file):
                    file_name = self.sample + "_CDS.fna"
                elif re.search(r"CDS\.gff", file):
                    self.logger.info("%s"%file)
                    file_name = self.sample + "_CDS.gff"
                elif re.search(r"rRNA\.fna", file):
                    file_name = self.sample + "_rRNA.fna"
                elif re.search(r"rRNA\.faa", file):
                    file_name = self.sample + "_rRNA.faa"
                elif re.search(r"rRNA\.gff", file):
                    file_name = self.sample + "_rRNA.gff"
                elif re.search(r"tRNA\.faa", file):
                    file_name = self.sample + "_tRNA.faa"
                elif re.search(r"tRNA\.fna", file):
                    file_name = self.sample + "_tRNA.fna"
                elif re.search(r"tRNA\.gff", file):
                    file_name = self.sample + "_tRNA.gff"
                out_name = os.path.join(last_path, file_name)
                os.system("cp {} {}".format(file_path, out_name))
            if os.path.exists(result_path):
                shutil.rmtree(result_path)
        else:
            link_dir(self.work_dir + "/merge_dir", last_path)
        if os.path.exists(last_path + "/" + self.sample + "_sample_stat.xls"):
            os.remove(last_path + "/" + self.sample + "_sample_stat.xls")
        os.link(self.stat.output_dir + "/" + self.sample + "_sample_stat.xls", last_path + "/" + self.sample + "_sample_stat.xls")
        if os.path.exists(last_path + "/" + self.sample + "_16S.fna"):
            os.remove(last_path + "/" + self.sample + "_16S.fna")
        if os.path.exists(self.stat.output_dir + "/" + self.sample + "_16S.fna"):
            if os.path.exists(last_path + "/" + self.sample + "_16S.fna"):
                os.remove(last_path + "/" + self.sample + "_16S.fna")
            os.link(self.stat.output_dir + "/" + self.sample + "_16S.fna", last_path + "/" + self.sample + "_16S.fna")
        if os.path.exists(result_path + "/" + self.sample + "_tRNA.faa"):
            os.remove(result_path + "/" + self.sample + "_tRNA.faa")
        if os.path.exists(result_path + "/" + self.sample + "_rRNA.faa"):
            os.remove(result_path + "/" + self.sample + "_rRNA.faa")
        statfile = self.work_dir + "/" + self.sample + "_base_stat.xls"
        statoutfile = last_path + "/" + self.sample + "_base_stat.xls"
        if os.path.exists(statfile):
            if os.path.exists(statoutfile):
                os.remove(statoutfile)
            os.link(statfile, statoutfile)
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行啦")
        super(GenePredictSingleModule, self).run()
        self.sample = self.get_params()
        if self.option("type") in ["majorbio"]:
            self.run_majorbio()
        elif self.option("type") in ["seq"]:
            self.run_seq()
        elif self.option("type") in ["gff"]:
            self.run_gff()
        if len(self.tools) == 1:
            self.on_rely(self.tools, self.predict_stat)
            for tool in self.tools:
                tool.run()
                gevent.sleep(0)
        elif len(self.tools) > 1:
            self.on_rely(self.tools, self.merge_file)
            for tool in self.tools:
                tool.run()
                gevent.sleep(0)

    def end(self):
        super(GenePredictSingleModule, self).end()

