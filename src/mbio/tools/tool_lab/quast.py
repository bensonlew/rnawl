# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# version 1.0
# last_modify: 2021.05.06

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
import json

class QuastAgent(Agent):
    """
    单个基因组QUAST基因组评估
    """
    def __init__(self, parent):
        super(QuastAgent, self).__init__(parent)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},# 序列文件夹
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta_dir").is_set:
            raise OptionError("必须输入序列文件夹")

    def set_resource(self):
        self._cpu = 8
        self._memory = '20G'

    def end(self):
        super(QuastAgent, self).end()

class QuastTool(Tool):
    def __init__(self, config):
        super(QuastTool, self).__init__(config)

    def run_quast(self):
        """
        运行phigaro
        :return:
        """
        cmd = "/bioinfo/tool_lab/Quast/quast-5.0.2/quast.py "
        if self.option("ref_fasta").is_set:
            cmd += "-R {} ".format(self.option("ref_fasta").prop["path"])
        cmd += "-o quast_ref "
        for fasta in os.listdir(self.option("fasta_dir").prop["path"]):
            if fasta.endswith(".fasta") or fasta.endswith(".fna") or fasta.endswith(".fa"):
                cmd += "{} ".format(self.option("fasta_dir").prop["path"] + "/" + fasta)
        self.logger.info(cmd)
        command = self.add_command("run_quast", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_quast运行完成！")
        else:
            self.set_error("run_quast运行出错!")

    def run_stat(self):
        """
        运行run_stat,主要统计信息
        :return:
        """
        if os.path.exists(self.output_dir + "/QUAST_Evaluation.xls"):
            os.remove(self.output_dir + "/QUAST_Evaluation.xls")
        with open(self.output_dir + "/QUAST_Evaluation.xls","w") as t, open(self.work_dir+"/quast_ref/report.tsv","r") as f:
            data = f.readlines()
            all_data_dict = {}
            for i in range(len(data[0].strip().split("\t")[1:])):
                all_data_dict[data[0].strip().split("\t")[i+1]] = {}
                for x in data[1:]:
                    if x.strip().split("\t")[0] == "# contigs":
                        all_data_dict[data[0].strip().split("\t")[i+1]]["contig"] = x.strip().split("\t")[i+1]
                    elif x.strip().split("\t")[0] == "Total length":
                        all_data_dict[data[0].strip().split("\t")[i+1]]["all_len"] = x.strip().split("\t")[i+1]
                    elif x.strip().split("\t")[0] == "Largest contig":
                        all_data_dict[data[0].strip().split("\t")[i+1]]["largest"] = x.strip().split("\t")[i+1]
                    elif x.strip().split("\t")[0] == "GC (%)":
                        all_data_dict[data[0].strip().split("\t")[i+1]]["gc"] = x.strip().split("\t")[i+1]
                    elif x.strip().split("\t")[0] == "N50":
                        all_data_dict[data[0].strip().split("\t")[i+1]]["n50"] = x.strip().split("\t")[i+1]
                    elif x.strip().split("\t")[0] == "N75":
                        all_data_dict[data[0].strip().split("\t")[i+1]]["n75"] = x.strip().split("\t")[i+1]
                    elif x.strip().split("\t")[0] == "L50":
                        all_data_dict[data[0].strip().split("\t")[i+1]]["l50"] = x.strip().split("\t")[i+1]
                    elif x.strip().split("\t")[0] == "L75":
                        all_data_dict[data[0].strip().split("\t")[i+1]]["l75"] = x.strip().split("\t")[i+1]
                    elif x.strip().split("\t")[0] == "# N's per 100 kbp":
                        all_data_dict[data[0].strip().split("\t")[i + 1]]["n_per"] = x.strip().split("\t")[i+1]
                    if self.option("ref_fasta").is_set:
                        if x.strip().split("\t")[0] == "Genome fraction (%)":
                            all_data_dict[data[0].strip().split("\t")[i+1]]["fraction"] = x.strip().split("\t")[i+1]
                        elif x.strip().split("\t")[0] == "Duplication ratio":
                            all_data_dict[data[0].strip().split("\t")[i+1]]["duplication"] = x.strip().split("\t")[i+1]
                        elif x.strip().split("\t")[0] == "Largest alignment":
                            all_data_dict[data[0].strip().split("\t")[i+1]]["alignment"] = x.strip().split("\t")[i+1]
                        elif x.strip().split("\t")[0] == "Total aligned length":
                            all_data_dict[data[0].strip().split("\t")[i+1]]["total_aligned"] = x.strip().split("\t")[i+1]
                        elif x.strip().split("\t")[0] == "# mismatches per 100 kbp":
                            all_data_dict[data[0].strip().split("\t")[i+1]]["mismatches"] = x.strip().split("\t")[i+1]
                        elif x.strip().split("\t")[0] == "# indels per 100 kbp":
                            all_data_dict[data[0].strip().split("\t")[i+1]]["indels"] = x.strip().split("\t")[i+1]
                        elif x.strip().split("\t")[0] == "Reference length":
                            all_data_dict[data[0].strip().split("\t")[i+1]]["ref_len"] = x.strip().split("\t")[i+1]
            if self.option("ref_fasta").is_set:
                t.write("Sample Name\tContigs\tTotal length（bp）\tLargest contig（bp）\tReference length\tGC Content (%)\tN50（bp）"
                        "\tL50\tN75（bp）\tL75\tN/100 kbp\tGenome fraction (%)\tDuplication ratio\tLargest alignment（bp）"
                        "\tTotal aligned length（bp）\tMismatches/100 kbp\tIndels/100 kbp\n")
            else:
                t.write("Sample Name\tContigs\tTotal length（bp）\tLargest contig（bp）\tGC Content (%)\tN50（bp）"
                    "\tL50\tN75（bp）\tL75\tN/100 kbp\n")
            for xx in all_data_dict:
                if self.option("ref_fasta").is_set:
                    t.write(xx + "\t" + all_data_dict[xx]["contig"] + "\t" + all_data_dict[xx]["all_len"]
                            + "\t" + all_data_dict[xx]["largest"] + "\t" + all_data_dict[xx]["ref_len"] + "\t" + all_data_dict[xx]["gc"]
                            + "\t" + all_data_dict[xx]["n50"] + "\t" + all_data_dict[xx]["l50"] + "\t" + all_data_dict[xx]["n75"]
                            + "\t" + all_data_dict[xx]["l75"] + "\t" + all_data_dict[xx]["n_per"] + "\t" + all_data_dict[xx]["fraction"]
                            + "\t" + all_data_dict[xx]["duplication"] + "\t" + all_data_dict[xx]["alignment"] + "\t" + all_data_dict[xx]["total_aligned"]
                            + "\t" + all_data_dict[xx]["mismatches"] + "\t" + all_data_dict[xx]["indels"] + "\n")
                else:
                    t.write(xx + "\t" + all_data_dict[xx]["contig"] + "\t" +
                            all_data_dict[xx]["all_len"] + "\t" + all_data_dict[xx]["largest"] + "\t" + all_data_dict[xx]["gc"]
                            + "\t" + all_data_dict[xx]["n50"] + "\t" + all_data_dict[xx]["l50"] + "\t" + all_data_dict[xx]["n75"]
                            + "\t" + all_data_dict[xx]["l75"] + "\t" + all_data_dict[xx]["n_per"] + "\n")

    def set_output(self):
        """
        设置结果文件目录
        """
        if os.path.exists(self.output_dir + "/cumulative_plot.pdf"):
            os.remove(self.output_dir + "/cumulative_plot.pdf")
        if self.option("ref_fasta").is_set:
            if os.path.exists(self.work_dir + "/quast_ref/aligned_stats/cumulative_plot.pdf"):
                os.link(self.work_dir + "/quast_ref/aligned_stats/cumulative_plot.pdf",
                        self.output_dir + "/cumulative_plot.pdf")
        else:
            if os.path.exists(self.work_dir + "/quast_ref/basic_stats/cumulative_plot.pdf"):
                os.link(self.work_dir + "/quast_ref/basic_stats/cumulative_plot.pdf",
                        self.output_dir + "/cumulative_plot.pdf")
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(QuastTool, self).run()
        self.run_quast()
        self.run_stat()
        self.set_output()
        self.end()