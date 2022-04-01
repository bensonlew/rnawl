# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20171218
"""高通量数据拆分样本拆分及质控一键化工作流"""


from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import json
import shutil
import re
import os


class SampleSplitQcWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        样本拆分及质控参数设置
        """
        self._sheet = wsheet_object
        super(SampleSplitQcWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "project_params", "type": "infile", "format": "datasplit.library_params"},  # 参数文件
            {'name': 'status_id', "type": 'string'},  # 分析记录表id
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.start_times = 0
        self.end_times = 0

    def check_options(self):
        if not self.option("project_params").is_set:
            raise OptionError("请设置参数文件")
        if not self.option("status_id"):
            raise OptionError("请设置状态表id")

    def get_params(self):
        """
        解析输入json文件，得到其中的参数信息
        """
        f = open(self.option("project_params").prop["path"], "rb")
        try:
            self.json_dict = json.loads(f.read())
        except:
            raise OptionError("json格式不正确")
        return self.json_dict

    def run_meta_split_qc(self):
        """
        运行多样性的质控及二次拆分的module
        """
        metas = []
        self.start_times += 1
        self.meta_start_times, self.meta_end_times = 0, 0
        meta_info = self.json_dict["meta"]
        for i in range(len(meta_info)):
            step_info = meta_info[i]
            options = {}
            options["lib_path"] = step_info["lib_path"]
            options["barcode_info"] = step_info["barcode_info"]
            options["lib_insert_size"] = step_info["lib_insert_size"]
            options["lib_specimen_id"] = step_info["lib_specimen_id"]
            types = step_info.keys()
            if "fq_type" in types:
                options["fq_type"] = step_info["fq_type"]
            if "leading" in types:
                options["leading"] = step_info["leading"]
            if "tailing" in types:
                options["tailing"] = step_info["tailing"]
            if "sliding_window" in types:
                options["sliding_window"] = step_info["sliding_window"]
            if "minlen" in types:
                options["minlen"] = step_info["minlen"]
            if "valid_len" in types:
                options["valid_len"] = step_info["valid_len"]
            if "min_lenth" in types:
                options["min_lenth"] = step_info["min_lenth"]
            if "max_lenth" in types:
                options["max_lenth"] = step_info["max_lenth"]
            if "mismatch_rate" in types:
                options["mismatch_rate"] = step_info["mismatch_rate"]
            if "pred" in types:
                options["pred"] = step_info["pred"]
            if "thread" in types:
                options["thread"] = step_info["thread"]
            if "min_len" in types:
                options["min_len"] = step_info["min_len"]
            if "split_type" in types:
                options["split_type"] = step_info["split_type"]
            self.meta_qc = self.add_module("datasplit.meta_qc")
            self.meta_qc.set_options(options)
            self.meta_qc.on('end', self.set_output, "meta_qc_{}".format(i))
            metas.append(self.meta_qc)
            self.meta_start_times += 1
        for m in metas:
            m.run()

    def run_metagenomic_qc(self):
        """
        运行宏基因组质控的module
        """
        metagenomics = []
        self.start_times += 1
        self.metagenomic_start_times, self.metagenomic_end_times = 0, 0
        metagenomic_info = self.json_dict["meta_genomic"]
        for i in range(len(metagenomic_info)):
            metagenomic_step = metagenomic_info[i]
            options = {}
            options["sample_path"] = metagenomic_step["sample_path"]
            if "quality_q" in metagenomic_step.keys():
                options["quality_q"] = metagenomic_step["quality_q"]
            if "length_q" in metagenomic_step.keys():
                options["length_q"] = metagenomic_step["length_q"]
            self.metagenomic_qc = self.add_module("datasplit.metagenomic_qc")
            self.metagenomic_qc.set_options(options)
            self.metagenomic_qc.on('end', self.set_output, "metagenomic_qc_{}".format(i))
            metagenomics.append(self.metagenomic_qc)
            self.metagenomic_start_times += 1
        for m in metagenomics:
            m.run()

    def run_microbial_genome_qc(self):
        """
        运行微生物基因组质控的module
        """
        microbial_genomes = []
        self.start_times += 1
        self.microbial_start_times, self.microbial_end_times = 0, 0
        microbial_info = self.json_dict["microbial_genome"]
        for i in range(len(microbial_info)):
            microbial_step = microbial_info[i]
            options = {}
            options["sample_path"] = microbial_step["sample_path"]
            options["sample_info"] = microbial_step["sample_info"]
            types = microbial_step.keys()
            if "readl" in types:
                options["readl"] = microbial_step["readl"]
            if "leading" in types:
                options["leading"] = microbial_step["leading"]
            if "tailing" in types:
                options["tailing"] = microbial_step["tailing"]
            if "sliding_window" in types:
                options["sliding_window"] = microbial_step["sliding_window"]
            if "minlen" in types:
                options["minlen"] = microbial_step["minlen"]
            if "seqprep_quality" in types:
                options["seqprep_quality"] = microbial_step["seqprep_quality"]
            if "seqprep_length" in types:
                options["seqprep_length"] = microbial_step["seqprep_length"]
            if "adapter_a" in types:
                options["adapter_a"] = microbial_step["adapter_a"]
            if "adapter_b" in types:
                options["adapter_b"] = microbial_step["adapter_b"]
            if "qual_type" in types:
                options["qual_type"] = microbial_step["qual_type"]
            self.microbial_qc = self.add_module("datasplit.microbial_genome_qc")
            self.microbial_qc.set_options(options)
            self.microbial_qc.on('end', self.set_output, "microbial_qc_{}".format(i))
            microbial_genomes.append(self.microbial_qc)
            self.microbial_start_times += 1
        for m in microbial_genomes:
            m.run()

    def run_mirna_qc(self):
        """
        运行miRNA质控的module
        """
        mirnas = []
        self.start_times += 1
        self.mirna_start_times, self.mirna_end_times = 0, 0
        mirna_info = self.json_dict["mirna"]
        for i in range(len(mirna_info)):
            mirna_step = mirna_info[i]
            options = {}
            options["list_file"] = mirna_step["list_file"]
            types = mirna_step.keys()
            if "length" in types:
                options["length"] = mirna_step["length"]
            if "adapter" in types:
                options["adapter"] = mirna_step["adapter"]
            if "phred_score" in types:
                options["phred_score"] = mirna_step["phred_score"]
            if "minlen" in types:
                options["minlen"] = mirna_step["minlen"]
            if "max_length" in types:
                options["max_length"] = mirna_step["max_length"]
            if "cut_left" in types:
                options["cut_left"] = mirna_step["cut_left"]
            self.mirna_qc = self.add_module("datasplit.mirna_qc")
            self.mirna_qc.set_options(options)
            self.mirna_qc.on('end', self.set_output, "mirna_qc_{}".format(i))
            mirnas.append(self.mirna_qc)
            self.mirna_start_times += 1
        for m in mirnas:
            m.run()

    def run_dna_split(self):
        """
        运行dna二次拆分的module
        """
        dna_splits = []
        self.start_times += 1
        self.split_start_times, self.split_end_times = 0, 0
        dna_info = self.json_dict["dna"]
        for i in range(len(dna_info)):
            split_step = dna_info[i]
            options = {}
            options["lib_path"] = split_step["lib_path"]
            options["library_info"] = split_step["lib_info"]
            if "ziplevel" in split_step.keys():
                options["ziplevel"] = split_step["ziplevel"]
            if "combinatorial" in split_step.keys():
                options["combinatorial"] = split_step["combinatorial"]
            self.dna_split = self.add_module("datasplit.dna_split_by_barcode")
            self.dna_split.set_options(options)
            self.dna_split.on("end", self.set_output, "dna_split_{}".format(i))
            dna_splits.append(self.dna_split)
            self.split_start_times += 1
        for m in dna_splits:
            m.run()

    def run_dna_rawdata_stat(self):
        """
        对dna拆分后的原始样本数据进行统计，运行原始数据统计的module
        """
        self.logger.info("开始进行dna样本原始数据的统计")
        list_file = os.path.join(self.work_dir, "dna_rawdata_list.txt")
        options = {
            "list_file": list_file
        }
        self.dna_rawdata_stat = self.add_module("datasplit.raw_data_stat")
        self.dna_rawdata_stat.set_options(options)
        self.dna_rawdata_stat.on("end", self.set_output, "dna_rawdata_stat")
        self.dna_rawdata_stat.run()

    def run_dna_qc(self):
        """
        运行dna质控的module:fastp
        """
        dnas = []
        self.start_times += 1
        self.dna_start_times, self.dna_end_times = 0, 0
        dna_info = self.json_dict["dna"]
        raw = open(self.work_dir + "/dna_rawdata_list.txt", "w")
        for i in range(len(dna_info)):
            dna_step = dna_info[i]
            options = {}
            types = dna_step.keys()
            if "qualified_quality_phred" in types:
                options["qualified_quality_phred"] = dna_step["qualified_quality_phred"]
            if "length_required" in types:
                options["length_required"] = dna_step["length_required"]
            if "cut_by_quality5" in types:
                options["cut_by_quality5"] = dna_step["cut_by_quality5"]
            if "cut_by_quality3" in types:
                options["cut_by_quality3"] = dna_step["cut_by_quality3"]
            if "cut_mean_quality" in types:
                options["cut_mean_quality"] = dna_step["cut_mean_quality"]
            if "n_base_limit" in types:
                options["n_base_limit"] = dna_step["n_base_limit"]
            if "compression" in types:
                options["compression"] = dna_step["compression"]
            if "thread" in types:
                options["thread"] = dna_step["thread"]
            sample_path = os.path.join(self.work_dir,  "dna_qc_list_{}.txt".format(i))
            fs = open(sample_path, "w")
            for lib in dna_step["samples"].keys():
                for s in dna_step["samples"][lib]:
                    for f in os.listdir(self.output_dir + "/dna_split"):
                        if re.match(r".*{}.*{}_R1.fastq.gz".format(lib, s), f):
                            f1 = self.output_dir + "/dna_split/" + lib + ":" + s + "_R1_raw.fastq.gz"
                            os.rename(self.output_dir + "/dna_split/" + f, f1)
                            fs.write(f1 + "\t" + lib + ":" + s + "\t" + "l\n")
                            raw.write(f1 + "\t" + lib + ":" + s + "\t" + "l\n")
                        elif re.match(r".*{}.*{}_R2.fastq.gz".format(lib, s), f):
                            f2 = self.output_dir + "/dna_split/" + lib + ":" + s + "_R2_raw.fastq.gz"
                            os.rename(self.output_dir + "/dna_split/" + f, f2)
                            fs.write(f2 + "\t" + lib + ":" + s + "\t" + "r\n")
                            raw.write(f2 + "\t" + lib + ":" + s + "\t" + "r\n")
            fs.close()
            options["sample_path"] = sample_path
            self.dna_qc = self.add_module("datasplit.fastp")
            self.dna_qc.set_options(options)
            self.dna_qc.on('end', self.set_output, "dna_qc_{}".format(i))
            dnas.append(self.dna_qc)
            self.dna_start_times += 1
        raw.close()
        # self.run_dna_rawdata_stat()  # 直接用fastp的json进行导表
        for m in dnas:
            m.run()

    # def run_rna_qc(self):
    #     """
    #     运行常规rna质控的module：fastp
    #     """
    #     rnas = []
    #     self.start_times += 1
    #     self.rna_start_times, self.rna_end_times = 0, 0
    #     rna_info = self.json_dict["rna"]
    #     for i in range(len(rna_info)):
    #         rna_step = rna_info[i]
    #         options = {}
    #         options["sample_path"] = rna_step["sample_path"]
    #         types = rna_step.keys()
    #         if "qualified_quality_phred" in types:
    #             options["qualified_quality_phred"] = rna_step["qualified_quality_phred"]
    #         if "length_required" in types:
    #             options["length_required"] = rna_step["length_required"]
    #         if "cut_by_quality5" in types:
    #             options["cut_by_quality5"] = rna_step["cut_by_quality5"]
    #         if "cut_by_quality3" in types:
    #             options["cut_by_quality3"] = rna_step["cut_by_quality3"]
    #         if "cut_mean_quality" in types:
    #             options["cut_mean_quality"] = rna_step["cut_mean_quality"]
    #         if "n_base_limit" in types:
    #             options["n_base_limit"] = rna_step["n_base_limit"]
    #         if "compression" in types:
    #             options["compression"] = rna_step["compression"]
    #         if "thread" in types:
    #             options["thread"] = rna_step["thread"]
    #         self.rna_qc = self.add_module("datasplit.fastp")
    #         self.rna_qc.set_options(options)
    #         self.step.add_steps("rna_qc_{}".format(i))
    #         step = getattr(self.step, 'rna_qc_{}'.format(i))
    #         step.start()
    #         self.step.update()
    #         self.rna_qc.on('end', self.finish_update, 'rna_qc_{}'.format(i))
    #         self.rna_qc.on('end', self.set_output, "rna_qc_{}".format(i))
    #         rnas.append(self.rna_qc)
    #         self.rna_start_times += 1
    #     for m in rnas:
    #         m.run()

    def run_rna_qc(self):  # 修改常规RNA质控流程为module: rna_qc 20181107
        """
        运行常规rna质控的module：rna_qc
        """
        rnas = []
        self.start_times += 1
        self.rna_start_times, self.rna_end_times = 0, 0
        rna_info = self.json_dict["rna"]
        for i in range(len(rna_info)):
            rna_step = rna_info[i]
            options = {}
            options["sample_path"] = rna_step["sample_path"]
            types = rna_step.keys()
            if "quality" in types:
                options["quality"] = rna_step["quality"]
            if "length" in types:
                options["length"] = rna_step["length"]
            if "adapter_a" in types:
                options["adapter_a"] = rna_step["adapter_a"]
            if "adapter_b" in types:
                options["adapter_b"] = rna_step["adapter_b"]
            self.rna_qc = self.add_module("datasplit.rna_qc")
            self.rna_qc.set_options(options)
            self.rna_qc.on('end', self.set_output, "rna_qc_{}".format(i))
            rnas.append(self.rna_qc)
            self.rna_start_times += 1
        for m in rnas:
            m.run()

    def run_ncrna_qc(self):
        """
        运行ncRNA质控的module
        """
        ncrnas = []
        self.start_times += 1
        self.ncrna_start_times, self.ncrna_end_times = 0, 0
        ncrna_info = self.json_dict["ncrna"]
        for i in range(len(ncrna_info)):
            ncrna_step = ncrna_info[i]
            options = {}
            options["list_file"] = ncrna_step["list_file"]
            types = ncrna_step.keys()
            if "fq_type" in types:
                options["fq_type"] = ncrna_step["fq_type"]
            if "minlen" in types:
                options["minlen"] = ncrna_step["minlen"]
            if "low_quality_base" in types:
                options["low_quality_base"] = ncrna_step["low_quality_base"]
            if "adaptor" in types:
                options["adaptor"] = ncrna_step["adaptor"]
            if "r_adaptor" in types:
                options["r_adaptor"] = ncrna_step["r_adaptor"]
            if "l_adaptor" in types:
                options["l_adaptor"] = ncrna_step["l_adaptor"]
            if "cut_left" in types:
                options["cut_left"] = ncrna_step["cut_left"]
            self.ncrna_qc = self.add_module("datasplit.ncrna_qc")
            self.ncrna_qc.set_options(options)
            self.ncrna_qc.on('end', self.set_output, "ncrna_qc_{}".format(i))
            ncrnas.append(self.ncrna_qc)
            self.ncrna_start_times += 1
        for m in ncrnas:
            m.run()

    def run_qc_stat(self, project_type):
        """
        运行质控后统计的module，统计质控后的结果
        """
        self.logger.info("开始进行{}质控后的统计".format(project_type))
        list_file = os.path.join(self.work_dir, project_type + "_stat_list.txt")
        fs = open(list_file, "w")
        qc_path = self.output_dir + "/" + project_type + "_qc"
        for f in os.listdir(qc_path):
            if project_type in ["ncrna", "mirna"]:
                if re.match(r".+xls", f):
                    continue
                m = re.match(r"(.+)_clean.fasta.gz", f)  # mirna样本
                n = re.match(r"(.+)_R[1,2]_clean.fastq.gz", f)  # ncrna样本
                if m:
                    s = m.group(1)
                elif n:
                    s = n.group(1)
                else:
                    raise Exception("没找到匹配的样本，请检查质控后的结果文件命名里是否有样本")
                if re.match(r".+R1.+", f):
                    fs.write(qc_path + "/" + f + "\t" + s + "\t" + "l\n")
                elif re.match(r".+R2.+", f):
                    fs.write(qc_path + "/" + f + "\t" + s + "\t" + "r\n")
                else:
                    fs.write(qc_path + "/" + f + "\t" + s + "\n")
            if project_type in ["rna", "dna", "microbial_genome", "meta", "meta_genomic"]:
                try:
                    s = f.split(".")[0]
                    if re.match(r".+clean.1.+", f):
                        fs.write(qc_path + "/" + f + "\t" + s + "\t" + "l\n")
                    elif re.match(r".+clean.2.+", f):
                        fs.write(qc_path + "/" + f + "\t" + s + "\t" + "r\n")
                    else:
                        if project_type == "microbial_genome":
                            s = s + "_unpaired"
                        fs.write(qc_path + "/" + f + "\t" + s + "\n")
                except:
                    raise Exception("dna/meta/microbial_genome没找到匹配的样本，请检查质控后的结果文件命名里是否有样本")
        fs.close()
        options = {
            "list_file": list_file,
            "project_type": project_type
        }
        self.qc_stat = self.add_module("datasplit.clean_data_stat")
        self.qc_stat.set_options(options)
        self.qc_stat.on("end", self.set_output, project_type + "_stat")
        self.qc_stat.run()

    def run_meta_qc_r1_stat(self):
        """
        运行质控后统计的module，统计meta只有R1端质控后的结果
        """
        self.logger.info("开始进行r1_meta质控后的统计")
        self.start_times += 1
        list_file = os.path.join(self.work_dir, "meta_stat_list_r1.txt")
        fs = open(list_file, "w")
        qc_path = self.output_dir + "/" + "meta_qc_r1"
        for f in os.listdir(qc_path):
            try:
                s = f.split(".fq.gz")[0]
                fs.write(qc_path + "/" + f + "\t" + s + "\n")
            except:
                raise Exception("没找到匹配的样本，请检查质控后的结果文件命名里是否有样本")
        fs.close()
        options = {
            "list_file": list_file,
            "project_type": "meta"
        }
        self.qc_stat = self.add_module("datasplit.clean_data_stat")
        self.qc_stat.set_options(options)
        self.qc_stat.on("end", self.set_output, "meta_stat_r1")
        self.qc_stat.run()

    def set_output(self, event):
        obj = event["bind_object"]
        self.logger.info("set output {}".format(event["data"]))
        if re.match(r"meta_qc_*", event["data"]):
            self.meta_end_times += 1
            for f in os.listdir(obj.output_dir):
                f1 = obj.output_dir + "/" + f + "/fastq_extract"
                self.links(f1, self.output_dir + "/meta_qc")
                f2 = obj.output_dir + "/" + f + "/R1_fastq_extract"
                if os.path.exists(f2):
                    self.links(f2, self.output_dir + "/meta_qc_r1")
            if self.meta_start_times == self.meta_end_times:
                self.run_qc_stat("meta")
                if os.path.exists(self.output_dir + "/meta_qc_r1"):
                    self.run_meta_qc_r1_stat()
        if re.match(r"metagenomic_qc_*", event["data"]):
            self.metagenomic_end_times += 1
            if not os.path.exists(self.output_dir + "/meta_genomic_qc"):
                os.mkdir(self.output_dir + "/meta_genomic_qc")
            for f in os.listdir(obj.output_dir):
                if f.endswith(".gz"):
                    if os.path.exists(self.output_dir + "/meta_genomic_qc/" + f):
                        os.remove(self.output_dir + "/meta_genomic_qc/" + f)
                    os.link(obj.output_dir + "/" + f, self.output_dir + "/meta_genomic_qc/" + f)
            if self.metagenomic_start_times == self.metagenomic_end_times:
                self.run_qc_stat("meta_genomic")
        if re.match(r"microbial_qc_*", event["data"]):
            self.microbial_end_times += 1
            for f in os.listdir(obj.output_dir):
                self.links(obj.output_dir, self.output_dir + "/microbial_genome_qc")
            if self.microbial_start_times == self.microbial_end_times:
                self.run_qc_stat("microbial_genome")
        if re.match(r"rna_qc_*", event["data"]):
            self.rna_end_times += 1
            self.links(obj.output_dir, os.path.join(self.output_dir, "rna_qc"))
            if self.rna_start_times == self.rna_end_times:
                self.run_qc_stat("rna")
        if re.match(r"mirna_qc_*", event["data"]):
            self.mirna_end_times += 1
            self.links(os.path.join(obj.output_dir, "fasta"), self.output_dir + "/mirna_qc")
            self.links(os.path.join(obj.output_dir, "fastq_stat"), self.output_dir + "/mirna_stat")
            if self.mirna_start_times == self.mirna_end_times:
                self.run_qc_stat("mirna")
        if re.match(r"ncrna_qc_*", event["data"]):
            self.ncrna_end_times += 1
            f2 = os.path.join(self.output_dir, "ncrna_qc")
            self.links(obj.output_dir, f2)
            if self.ncrna_start_times == self.ncrna_end_times:
                self.run_qc_stat("ncrna")
        if re.match(r"dna_split.*", event["data"]):
            self.split_end_times += 1
            self.links(obj.output_dir, os.path.join(self.output_dir, "dna_split"))
            if self.split_start_times == self.split_end_times:
                self.run_dna_qc()
                # self.run_api(project_type="dna_split", dir_path=os.path.join(self.output_dir, "dna_split"))
        if re.match(r"dna_qc.*", event["data"]):
            self.dna_end_times += 1
            self.links(obj.output_dir + "/fastq", os.path.join(self.output_dir, "dna_qc"))
            self.links(obj.output_dir + "/json_dir", os.path.join(self.output_dir, "dna_data_stat"))
            if self.dna_start_times == self.dna_end_times:
                # self.run_qc_stat("dna")  # 直接用fastp的json进行导表
                self.run_api(project_type="dna_fastp", dir_path=os.path.join(self.output_dir, "dna_data_stat"))
        # if re.match(r"dna_rawdata_stat", event["data"]):
        #     self.links(obj.output_dir, os.path.join(self.output_dir, "dna_rawdata_stat"))
        #     self.run_api(project_type="dna_rawdata_stat", dir_path=os.path.join(self.output_dir, "dna_rawdata_stat"))
        if re.match(r"rna_stat", event["data"]):
            self.links(obj.output_dir, os.path.join(self.output_dir, "rna_stat"))
            self.run_api(project_type="rna_stat", dir_path=os.path.join(self.output_dir, "rna_stat"))
        if re.match(r"mirna_stat", event["data"]):
            self.links(obj.output_dir, os.path.join(self.output_dir, "mirna_stat"))
            self.run_api(project_type="mirna_stat", dir_path=os.path.join(self.output_dir, "mirna_stat"))
        if re.match(r"ncrna_stat", event["data"]):
            self.links(obj.output_dir, os.path.join(self.output_dir, "ncrna_stat"))
            self.run_api(project_type="ncrna_stat", dir_path=os.path.join(self.output_dir, "ncrna_stat"))
        if re.match(r"dna_stat", event["data"]):
            self.links(obj.output_dir, os.path.join(self.output_dir, "dna_stat"))
            self.run_api(project_type="dna_stat", dir_path=os.path.join(self.output_dir, "dna_stat"))
        if event["data"] == "meta_stat":
            self.links(obj.output_dir, os.path.join(self.output_dir, "meta_stat"))
            self.run_api(project_type="meta_stat", dir_path=os.path.join(self.output_dir, "meta_stat"))
        if event["data"] == "meta_stat_r1":
            self.links(obj.output_dir, os.path.join(self.output_dir, "meta_stat_r1"))
            self.run_api(project_type="meta_stat_r1", dir_path=os.path.join(self.output_dir, "meta_stat_r1"))
        if re.match(r"meta_genomic_stat", event["data"]):
            self.links(obj.output_dir, os.path.join(self.output_dir, "meta_genomic_stat"))
            self.run_api(project_type="meta_genomic_stat", dir_path=os.path.join(self.output_dir, "meta_genomic_stat"))
        if re.match(r"microbial_genome_stat", event["data"]):
            self.links(obj.output_dir, os.path.join(self.output_dir, "microbial_genome_stat"))
            self.run_api(project_type="microbial_genome", dir_path=os.path.join(self.output_dir, "microbial_genome_stat"))

    def links(self, old, new):
        if not os.path.exists(new):
            os.mkdir(new)
        for f in os.listdir(old):
            if re.match(r"list.txt", f):
                with open(new + "/list.txt", "w+") as w, open(old + "/list.txt", "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        w.write(line)
                continue
            if os.path.isfile(os.path.join(old, f)):
                f_ = os.path.join(new, f)
                if os.path.exists(f_):
                    os.remove(f_)
                os.link(os.path.join(old, f), f_)
            else:
                if not os.path.exists(os.path.join(new, f)):
                        os.mkdir(os.path.join(new, f))
                for f1 in os.listdir(os.path.join(old, f)):
                    f_new = os.path.join(os.path.join(new, f), f1)
                    if os.path.exists(f_new):
                        os.remove(f_new)
                    os.link(os.path.join(os.path.join(old, f), f1), f_new)

    def run_api(self, project_type, dir_path):
        """
        对二次拆分及质控结果进行导表
        """
        self.logger.info("开始进行{}导表".format(project_type))
        self.end_times += 1
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        datasplit = self.api.datasplit_new
        if project_type == "dna_fastp":
            self.end_times += 1
            datasplit.add_fastp_stat(status_id=self.option("status_id"), json_dir=dir_path)
        # if project_type == "dna_rawdata_stat":
        #     datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="raw_sample", project_type="dna", dir_path=dir_path)
        if project_type == "ncrna_stat":
            datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="clean_sample", project_type="ncrna", dir_path=dir_path)
        if project_type == "rna_stat":
            datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="clean_sample", project_type="rna", dir_path=dir_path)
        if project_type == "mirna_stat":
            datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="clean_sample", project_type="mirna", dir_path=dir_path)
        # if project_type == "dna_stat":
        #     datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="clean_sample", project_type="dna", dir_path=dir_path)
        if project_type == "meta_stat":
            datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="clean_sample", project_type="meta", dir_path=dir_path)
        if project_type == "meta_stat_r1":
            datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="clean_sample", project_type="meta_r1", dir_path=dir_path)
        if project_type == "meta_genomic_stat":
            datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="clean_sample", project_type="meta_genomic", dir_path=dir_path)
        if project_type == "microbial_genome":
            datasplit.add_raw_stat(status_id=self.option("status_id"), data_type="clean_sample", project_type="microbial_genome", dir_path=dir_path)
        self.logger.info(self.start_times)
        self.logger.info(self.end_times)
        if self.start_times == self.end_times:
            if self._sheet.output.startswith("s3:"):
                s3_upload_dir = self._sheet.output
            else:
                s3_upload_dir = None
            datasplit.update_specimen_clean_path(status_id=self.option("status_id"), dir_path=self.output_dir, s3_upload_dir=s3_upload_dir)
            self.end()

    def run(self):
        self.get_params()
        if "meta" in self.json_dict.keys():
            self.run_meta_split_qc()
        if "meta_genomic" in self.json_dict.keys():
            self.run_metagenomic_qc()
        if "microbial_genome" in self.json_dict.keys():
            self.run_microbial_genome_qc()
        if "rna" in self.json_dict.keys():
            self.run_rna_qc()
        if "mirna" in self.json_dict.keys():
            self.run_mirna_qc()
        if "ncrna" in self.json_dict.keys():
            self.run_ncrna_qc()
        if "dna" in self.json_dict.keys():
            self.run_dna_split()
        super(SampleSplitQcWorkflow, self).run()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SampleSplitQcWorkflow, self).end()
