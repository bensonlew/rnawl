# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# version 1.0
# last_modify: 2018.03.06

import os, re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna.trans_step import step_count  # 对fasta序列做长度分布图输入表


class SampleGeneStatAgent(Agent):
    """
    样品层面的编码预测结果统计和基因序列长度分布统计（完成图的染色体和质粒序列合并成一个文件统计）
    """

    def __init__(self, parent):
        super(SampleGeneStatAgent, self).__init__(parent)
        options = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测结果
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "repeat_gff", "type": "infile", "format": "gene_structure.gff3"},  # repeat预测结果
            {"name": "seq", "type": "infile", "format": "sequence.fasta"},  # 染色体序列
            {"name": "gene_faa", "type": "infile", "format": "sequence.fasta"},
            {"name": "trna_fnn", "type": "infile", "format": "sequence.fasta"},
            {"name": "rrna_fnn", "type": "infile", "format": "sequence.fasta"},
            {"name": "repeat_dat", "type": "string", 'default': ""},
            {"name": "trna_struc", "type": "string", 'default': ""},
            {"name": "genome_plas", "type": "infile", "format": "sequence.fasta"},
            {"name": "gene_plas_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测结果
            {"name": "trna_plas_gff", "type": "infile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "rrna_plas_gff", "type": "infile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "repeat_plas_gff", "type": "infile", "format": "gene_structure.gff3"},  # repeat预测结果
            {"name": "seq_plas", "type": "infile", "format": "sequence.fasta"},  # 完成图的质粒序列
            {"name": "gene_plas_faa", "type": "infile", "format": "sequence.fasta"},
            {"name": "trna_plas_fnn", "type": "infile", "format": "sequence.fasta"},
            {"name": "rrna_plas_fnn", "type": "infile", "format": "sequence.fasta"},
            {"name": "repeat_plas_dat", "type": "string", 'default': ""},
            {"name": "trna_plas_struc", "type": "string", 'default': ""},
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "sample_gene_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 样品编码基因预测gff文件
            {"name": "sample_trna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "sample_rrna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "sample_repeat_gff", "type": "outfile", "format": "gene_structure.gff3"},  # TRE预测结果
            {"name": "sample_gene_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample_gene_faa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            {"name": "length_distribu", "type": "outfile", "format": "paternity_test.tab"},  # 样品基因序列的长度分布文件

        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("genome").is_set:
            raise OptionError("必须提供基因组文件", code="33300801")
        if not self.option("seq").is_set:
            raise OptionError("必须提供预测得到的基因序列！", code="33300802")
        if self.option("genome_plas").is_set and not self.option("seq_plas").is_set:
            raise OptionError("必须提供genome_plas预测得到的基因序列！", code="33300803")
        if not self.option("genome_plas").is_set and self.option("seq_plas").is_set:
            raise OptionError("必须提供预测基因序列的genome_plas！", code="33300804")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(SampleGeneStatAgent, self).end()


class SampleGeneStatTool(Tool):
    def __init__(self, config):
        super(SampleGeneStatTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.perl_script2 = self.config.PACKAGE_DIR + "/bacgenome/"

    def run_gene_info_stat(self):
        if self.option("genome").is_set and self.option("genome_plas").is_set:
            all_fasta = self.work_dir + "/genome.fna"
            os.system('cat ' + self.option("genome").prop['path'] + " " + self.option("genome_plas").prop[
                'path'] + ' >' + all_fasta)
        else:
            all_fasta = self.option("genome").prop['path']
        if self.option("seq").is_set and self.option("seq_plas").is_set:
            all_gene_fnn = self.work_dir + "/sample_all.fnn"
            os.system('cat ' + self.option("seq").prop['path'] + " " + self.option("seq_plas").prop[
                'path'] + ' >' + all_gene_fnn)
        elif self.option("seq").is_set and not self.option("seq_plas").is_set:
            all_gene_fnn = self.option("seq").prop['path']
        self.logger.info(step_count(all_gene_fnn, self.work_dir + "/fnn_stat.xls", 11, 100,
                                    self.output_dir + "/length_distribute.txt"))
        path = self.output_dir + "/" + "CDS_predict"
        if not os.path.exists(path):
            os.makedirs(path)
        step_count(all_gene_fnn, self.work_dir + "/fnn_stat.xls", 11, 100, self.work_dir + "/length_distribute.txt")
        cmd = '{} {}dnabac_sample_stat.pl {} {} {} {}'.format(self.perl_path, self.perl_script, all_fasta, all_gene_fnn,
                                                              self.option("sample"), path)
        command = self.add_command("sample_gene_info_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("样品编码基因预测统计完成")
        else:
            self.set_error("样品编码基因预测统计出错!", code="33300801")

    def file_merge_split(self, file1, file2, soft_type, file_type):
        path = ""
        if soft_type == "CDS":
            path = self.output_dir + "/" + "CDS_predict"
        elif soft_type == "TRF":
            path = self.output_dir + "/" + "repeats"
        else:
            path = self.output_dir + "/" + soft_type
        if not os.path.exists(path):
            os.makedirs(path)
        if file_type == "gff":
            f1 = pd.read_table(file1, sep='\t', header=0)
            if file2 != "":
                f2 = pd.read_table(file2, sep='\t', header=0)
                summary = pd.concat([f1, f2])
                f1 = summary
            if len(f1) < 1:
                f1.to_csv(path + "/" + self.option("sample") + "_" + soft_type + ".gff", sep="\t",
                          index=None)
                return
            if soft_type == "TRF":
                f1["location"] = f1["Sequence Name"]
            else:
                if len(f1["Sequence id"]) > 0:
                    f1["location"] = f1["Sequence id"].str.split("_", expand=True)[0]
            if re.match("^[Ss]caffold[0-9]+$", list(f1["location"])[0]) or soft_type == "rRNA" or soft_type == "tRNA":
                del f1["location"]
                f1.to_csv(path + "/" + self.option("sample") + "_" + soft_type + ".gff", sep="\t",
                          index=None)
            else:
                for loc in list(set(f1['location'])):
                    ta = f1[f1["location"] == loc]
                    location = loc
                    del ta["location"]
                    ta.to_csv(path + "/" + self.option("sample") + "_" + location + "_" + soft_type + ".gff",
                              sep="\t", index=None)
                del f1["location"]
                f1.to_csv(path + "/" + self.option("sample") + "_whole_genome_" + soft_type + ".gff",
                          sep="\t", index=None)
        else:
            all_gene_seq = self.work_dir + "/" + soft_type + "." + file_type
            os.system('cat ' + file1 + " " + file2 + ' >' + all_gene_seq)
            if os.path.getsize(all_gene_seq):
                with open(all_gene_seq, "r") as f:
                    if re.match("^.*[Ss]caffold[0-9].*", f.readlines()[0]):
                        if os.path.exists(path + "/" + self.option("sample") + "_" + soft_type + "." + file_type):
                            os.remove(path + "/" + self.option("sample") + "_" + soft_type + "." + file_type)
                        os.link(all_gene_seq, path + "/" + self.option("sample") + "_" + soft_type + "." + file_type)
                    else:
                        if file_type == "faa":
                            if os.path.exists(path + "/" + self.option("sample") + "_whole_genome_CDS" + ".faa"):
                                os.remove(path + "/" + self.option("sample") + "_whole_genome_CDS" + ".faa")
                            os.link(all_gene_seq,
                                    path + "/" + self.option("sample") + "_whole_genome_CDS" + ".faa")
                        if file_type == "fnn" and soft_type == "CDS":
                            if os.path.exists(path + "/" + self.option("sample") + "_whole_genome_CDS" + ".fnn"):
                                os.remove(path + "/" + self.option("sample") + "_whole_genome_CDS" + ".fnn")
                            os.link(all_gene_seq,
                                    path + "/" + self.option("sample") + "_whole_genome_CDS" + ".fnn")
                        if file_type == "dat" and soft_type == "TRF":
                            if os.path.exists(path + "/" + self.option("sample") + "_whole_genome_TRF" + ".dat"):
                                os.remove(path + "/" + self.option("sample") + "_whole_genome_TRF" + ".dat")
                            os.link(all_gene_seq, path + "/" + self.option("sample") + "_whole_genome_TRF" + ".dat")
                        elif file_type != "struc" and soft_type not in ["CDS", "TRF"]:
                            if os.path.exists(path + "/" + self.option("sample") + "_" + soft_type + "." + file_type):
                                os.remove(path + "/" + self.option("sample") + "_" + soft_type + "." + file_type)
                            os.link(all_gene_seq,
                                    path + "/" + self.option("sample") + "_" + soft_type + "." + file_type)
                        if soft_type in ["CDS", "TRF"] or file_type in ["struc"]:
                            cmd = '{} {}file_split_byloc.pl {} {} {}'.format(self.perl_path, self.perl_script2,
                                                                             all_gene_seq,
                                                                             self.option("sample"), path)
                            command = self.add_command(soft_type.lower() + file_type.lower(), cmd).run()
                            self.wait(command)
                            if command.return_code == 0:
                                self.logger.info("file_split_byloc完成")
                            else:
                                self.set_error("file_split_byloc出错!", code="33300802")
            else:
                return

    def set_output(self):
        sample_stat = self.output_dir + "/CDS_predict/" + self.option("sample") + "_CDS_statistics.xls"
        if os.path.exists(sample_stat):
            self.option('gene_statistics', sample_stat)
        else:
            self.option('gene_statistics',
                        self.output_dir + "/CDS_predict/" + self.option("sample") + "_whole_genome_CDS_statistics.xls")
        self.option('length_distribu', self.work_dir + "/length_distribute.txt")
        self.logger.info(self.option('gene_statistics').prop['path'])
        self.logger.info(self.option('length_distribu').prop['path'])
        cds_gff = self.output_dir + "/CDS_predict/" + self.option("sample") + "_whole_genome_CDS.gff"
        if os.path.exists(cds_gff):
            self.option('sample_gene_gff', cds_gff)
        else:
            self.option('sample_gene_gff', self.output_dir + "/" + "CDS_predict/" + self.option("sample") + "_CDS.gff")
        cds_faa = self.output_dir + "/CDS_predict/" + self.option("sample") + "_whole_genome_CDS.faa"
        if os.path.exists(cds_faa):
            self.option('sample_gene_faa', cds_faa)
        else:
            self.option('sample_gene_faa', self.output_dir + "/" + "CDS_predict/" + self.option("sample") + "_CDS.faa")
        cds_fnn = self.output_dir + "/CDS_predict/" + self.option("sample") + "_whole_genome_CDS.fnn"
        if os.path.exists(cds_fnn):
            self.option('sample_gene_fnn', cds_fnn)
        else:
            self.option('sample_gene_fnn', self.output_dir + "/CDS_predict/" + self.option("sample") + "_CDS.fnn")
        rrna_gff = self.output_dir + "/rRNA/" + self.option("sample") + "_rRNA.gff"
        if os.path.exists(rrna_gff):
            self.option('sample_rrna_gff', rrna_gff)
        trna_gff = self.output_dir + "/tRNA/" + self.option("sample") + "_tRNA.gff"
        if os.path.exists(trna_gff):
            self.option('sample_trna_gff', trna_gff)
        repeat_gff = self.output_dir + "/repeats/" + self.option("sample") + "_whole_genome_TRF.gff"
        if os.path.exists(repeat_gff):
            self.option('sample_repeat_gff', trna_gff)
        else:
            self.option('sample_repeat_gff', self.output_dir + "/repeats/" + self.option("sample") + "_TRF.gff")

    def run(self):
        super(SampleGeneStatTool, self).run()
        self.run_gene_info_stat()
        if self.option("gene_gff").is_set and self.option("gene_plas_gff").is_set:
            self.file_merge_split(self.option("gene_gff").prop['path'], self.option("gene_plas_gff").prop['path'],
                                  "CDS", "gff")
        elif self.option("gene_gff").is_set:
            self.file_merge_split(self.option("gene_gff").prop['path'], "", "CDS", "gff")
        if self.option("rrna_gff").is_set and self.option("rrna_plas_gff").is_set:
            self.file_merge_split(self.option("rrna_gff").prop['path'], self.option("rrna_plas_gff").prop['path'],
                                  "rRNA", "gff")
        elif self.option("rrna_gff").is_set:
            self.file_merge_split(self.option("rrna_gff").prop['path'], "", "rRNA", "gff")
        if self.option("trna_gff").is_set and self.option("trna_plas_gff").is_set:
            self.file_merge_split(self.option("trna_gff").prop['path'], self.option("trna_plas_gff").prop['path'],
                                  "tRNA", "gff")
        elif self.option("trna_gff").is_set:
            self.file_merge_split(self.option("trna_gff").prop['path'], "", "tRNA", "gff")
        if self.option("repeat_gff").is_set and self.option("repeat_plas_gff").is_set:
            self.file_merge_split(self.option("repeat_gff").prop['path'], self.option("repeat_plas_gff").prop['path'],
                                  "TRF", "gff")
        elif self.option("repeat_gff").is_set:
            self.file_merge_split(self.option("repeat_gff").prop['path'], "", "TRF", "gff")
        if self.option("seq").is_set and self.option("seq_plas").is_set:
            self.file_merge_split(self.option("seq").prop['path'], self.option("seq_plas").prop['path'], "CDS", "fnn")
        elif self.option("seq").is_set:
            self.file_merge_split(self.option("seq").prop['path'], "", "CDS", "fnn")
        if self.option("gene_faa").is_set and self.option("gene_plas_faa").is_set:
            self.file_merge_split(self.option("gene_faa").prop['path'], self.option("gene_plas_faa").prop['path'],
                                  "CDS", "faa")
        elif self.option("gene_faa").is_set:
            self.file_merge_split(self.option("gene_faa").prop['path'], "", "CDS", "faa")
        if self.option("rrna_fnn").is_set and self.option("rrna_plas_fnn").is_set:
            self.file_merge_split(self.option("rrna_fnn").prop['path'], self.option("rrna_plas_fnn").prop['path'],
                                  "rRNA", "fnn")
        elif self.option("rrna_fnn").is_set:
            self.file_merge_split(self.option("rrna_fnn").prop['path'], "", "rRNA", "fnn")
        if self.option("trna_fnn").is_set and self.option("trna_plas_fnn").is_set:
            self.file_merge_split(self.option("trna_fnn").prop['path'], self.option("trna_plas_fnn").prop['path'],
                                  "tRNA", "fnn")
        elif self.option("trna_fnn").is_set:
            self.file_merge_split(self.option("trna_fnn").prop['path'], "", "tRNA", "fnn")
        if self.option("repeat_dat") != "" and self.option("repeat_plas_dat") != "":
            self.file_merge_split(self.option("repeat_dat"), self.option("repeat_plas_dat"), "TRF", "dat")
        elif self.option("repeat_dat") != "":
            self.file_merge_split(self.option("repeat_dat"), "", "TRF", "dat")
        if self.option("trna_struc") != "" and self.option("trna_plas_struc") != "":
            self.logger.info(self.option("trna_struc"))
            self.logger.info(self.option("trna_plas_struc"))
            self.file_merge_split(self.option("trna_struc"), self.option("trna_plas_struc"), "tRNA", "struc")
        elif self.option("trna_struc") != "":
            self.file_merge_split(self.option("trna_struc"), "", "tRNA", "struc")
        self.set_output()
        self.end()
