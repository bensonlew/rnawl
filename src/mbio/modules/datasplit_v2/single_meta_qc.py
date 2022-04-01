# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190125

import os
import re
import json
import time
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class SingleMetaQcModule(Module):
    """
    多样性对单个文库进行二次拆分及质控
    """
    def __init__(self, work_id):
        super(SingleMetaQcModule, self).__init__(work_id)
        options = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 输入文件, 输入的单个文库文件夹，存放序列两个gz文件
            {"name": "barcode_info", "type": "infile", "format": "datasplit.meta_barcode_info"},  # 文库中样本barcode信息
            {"name": "primer_info", "type": "infile", "format": "datasplit.meta_primer_info"},  # 文库中样引物信息
            {"name": "sample_primer", "type": "infile", "format": "datasplit.path"},  # 文库中样本primer信息
            {"name": "lib_name", "type": "string"},  # 文库编号
            {"name": "lib_insert_size", "type": "int"},  # 文库插入片段长度
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE or SE
            {"name": "leading", "type": "string", "default": "0"},  # 切除首端碱基质量小于0的碱基或者N
            {"name": "tailing", "type": "string", "default": "20"},  # 切除末端碱基质量小于20的碱基或者N
            {"name": "sliding_window", "type": "string", "default": "50:20"}, # 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {"name": "minlen", "type": "string", "default": "50"},  # 最低reads长度
            {"name": "valid_len", "type": "int"},  # -l,长度过滤阈值
            {"name": "min_lenth", "type": "string", "default": "10"},  # -m,两个reads之间所需的最小重叠长度，以提供可靠的重叠
            {"name": "max_lenth", "type": "string", "default": "100"},  # -M,两个reads之间的最大重叠长度
            {"name": "mismatch_rate", "type": "string", "default": "0.2"},  # -x,错配和重叠长度允许的最大比率
            {"name": "pred", "type": "string", "default": "33"},  # -p,FASTQ文件中的碱基的质量值，Pred33/Pred64.
            {"name": "thread", "type": "string", "default": "6"},  # -t,线程数
            {"name": "min_len", "type": "int"},  # -m,最小长度
            {"name": "split_type", "type": "string", "default": "Auto"},  # 拆分样本序列类型 Pair or Single or Auto
            {"name": "its_primer", "type": "bool", "default": False},  # 是否是ITS引物，ITS引物质控前加SeqPrep
            # modifed by zengjing 20200102 SeqPrep、Trimmomatic的质控换成fastp
            {"name": "length_required", "type": "string", "default": "50"},  # -l,长度过滤参数，比此值短的读取将被丢弃，默认15
            {"name": "cut_right_mean_quality", "type": "string", "default": "20"},  # cut_right的平均质量要求,默认20
            {"name": "cut_right_window_size", "type": "string", "default": "50"},  # cut_right的窗口大小，默认4
            {"name": "cut_by_quality5", "type": "string", "default": "0"},  # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {"name": "cut_by_quality3", "type": "string", "default": "20"},  #  -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {"name": "lib_type", "type": "string", "default": "no_official"},  # 文库类型，官方建库和非官方建库
            # modifed by xueqinwen 20211011 增加允许引物错配数优化
            {'name': "mismatch", "type": "string"},  # -l，允许引物错配数
        ]
        self.add_option(options)
        self.fq_list = []  # 文库中的序列R1端/R2端
        self.valid_stat, self.fastp_json, self.flash_cmd, self.split_fq = "", "", "", "" # 文库质控统计
        self.raw_fastq, self.raw_fastq_end, self.fastq_extract_end = False, False, False
        if self.option("lib_type") != "official":
            self.raw_fastq = True

    def check_options(self):
        if not self.option("fq_dir").is_set:
            raise OptionError("请设置文库fastq文件夹")
        if not self.option("barcode_info").is_set:
            raise OptionError("请设置文库的样本信息表")
        if not self.option("primer_info").is_set:
            raise OptionError("请设置文库的引物信息表")
        if not self.option("sample_primer").is_set:
            raise OptionError("请设置样本的引物信息表")
        if not self.option("lib_insert_size"):
            raise OptionError("请设置文库插入片段长度")
        if not self.option("lib_name"):
            raise OptionError("请设置文库编号")
        if self.option("split_type") not in ["Pair", "Single", "Auto"]:
            raise OptionError("拆分类型只能是Pair或Single或自动拆分")
        if self.option("lib_type") not in ["official", "no_official"]:
            raise OptionError("文库类型只能是official或no_official")

    def check_primer_info(self):
        barcode_list = []
        with open(self.option("primer_info").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                seq = item[1] + "_" + item[3]
                if seq in barcode_list:
                    self.logger.info(self.option("lib_type"))
                    if self.option("lib_type") != "official":
                        raise OptionError("{}里的barcode：{}重复，请检查".format(self.option("primer_info").prop["path"], seq))
                barcode_list.append(seq)

    def check_barcode_info(self):
        barcode_name_list, barcode_seq_list = [], []
        with open(self.option("barcode_info").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                if self.option("lib_type") != "official":
                    if item[1] in barcode_name_list:
                        raise OptionError("{}里的barcode：{}重复，请检查".format(self.option("barcode_info").prop["path"], item[1]))
                    barcode_name_list.append(item[1])
                    seq = item[2] + "_" + item[3]
                    if seq in barcode_seq_list:
                        raise OptionError("{}里的barcode：{}重复，请检查".format(self.option("barcode_info").prop["path"], seq))
                    barcode_seq_list.append(seq)

    def run_getvaild_bybarcode(self):
        options = ({
            "fq1": self.fq_list[0],
            "fq2": self.fq_list[1],
            "barcode_info": self.option("barcode_info"),
            "lib_name": self.option("lib_name"),
        })
        self.getvaild_bybarcode = self.add_tool("datasplit_v2.getvaild_bybarcode")
        self.getvaild_bybarcode.set_options(options)
        self.getvaild_bybarcode.on("end", self.run_fastp)
        self.getvaild_bybarcode.on("end", self.run_get_raw_fastq, "valid")
        self.getvaild_bybarcode.on("end", self.set_output, "valid")
        self.getvaild_bybarcode.run()

    def run_get_raw_fastq(self, event):
        if self.option("lib_type") != "official":
            obj = event["bind_object"]
            fq1 = self.getvaild_bybarcode.option("out_fq1").prop["path"]
            fq2 = self.getvaild_bybarcode.option("out_fq2").prop["path"]
            for f in os.listdir(obj.work_dir):
                if f.endswith("seq2sam.stat"):
                    seq2sam = os.path.join(obj.work_dir, f)
            self.logger.info(seq2sam)
            is_end = False
            if os.path.getsize(fq1) == 0:
                is_end = True
            elif os.path.getsize(fq2) == 0:
                is_end = True
            if not is_end:
                options = {
                    # "fq1": self.fq_list[0],
                    # "fq2": self.fq_list[1],
                    "fq1": fq1,
                    "fq2": fq2,
                    "seq2sam": seq2sam,
                    "sample_primer": self.option("sample_primer"),
                }
                self.fastq_extract_raw = self.add_tool("datasplit_v2.fastq_extract_raw")
                self.fastq_extract_raw.set_options(options)
                self.fastq_extract_raw.on("end", self.set_output, "raw_fastq")
                self.fastq_extract_raw.run()

    def run_fastp(self):
        """
        fastp进行质控
        """
        if self.option("lib_type") == "official":
            fq1 = self.fq_list[0]
            fq2 = self.fq_list[1]
        else:
            fq1 = self.getvaild_bybarcode.option("out_fq1").prop["path"]
            fq2 = self.getvaild_bybarcode.option("out_fq2").prop["path"]
        if os.path.getsize(fq1) == 0:
            self.end()
        elif os.path.getsize(fq2) == 0:
            self.end()
        else:
            options = {
                "fq1": fq1,
                "fq2": fq2,
                "length_required": self.option("length_required"),
                "cut_right_mean_quality": self.option("cut_right_mean_quality"),
                "cut_right_window_size": self.option("cut_right_window_size"),
                "cut_by_quality5": self.option("cut_by_quality5"),
                "cut_by_quality3": self.option("cut_by_quality3")
            }
            self.fastp = self.add_tool("datasplit_v2.fastp_meta")
            self.fastp.set_options(options)
            if self.split_type != "Single":  # 不是单端拆分，则需要拼接
                if self.option("lib_insert_size") <= 380 and not self.option("its_primer"):
                    self.fastp.on("end", self.run_trim_length_flash)
                else:
                    self.fastp.on("end", self.run_flash)
            else:
                self.fastp.on("end", self.run_R1_split_by_barcode)
            self.fastp.on("end", self.set_output, "fastp")
            self.fastp.run()

    def run_trim_length_flash(self):  # 拼接前长度过滤
        self.tools, fq_list = [], []
        for f in os.listdir(self.fastp.output_dir):
            if f.endswith("trim.1.fq"):
                fq_list.insert(0, os.path.join(self.fastp.output_dir, f))
            if f.endswith("trim.2.fq"):
                fq_list.append(os.path.join(self.fastp.output_dir, f))
        for fq in fq_list:
            options = {
                "fq": fq,
                "valid_len": (str(self.option("valid_len")) if self.option("valid_len") else str(self.trim_length)),
                "lib_name": self.option("lib_name"),
            }
            self.trim_length_flash = self.add_tool("datasplit_v2.trim_length")  # 拼接前的长度过滤
            self.trim_length_flash.set_options(options)
            self.tools.append(self.trim_length_flash)
        self.on_rely(self.tools, self.run_flash)
        for tool in self.tools:
            tool.run()

    def run_flash(self):
        options = {
            "min_lenth": self.option("min_lenth"),
            "max_lenth": self.option("max_lenth"),
            "mismatch_rate": self.option("mismatch_rate"),
            "pred": self.option("pred"),
            "thread": self.option("thread"),
            "lib_name": self.option("lib_name")
        }
        if self.option("lib_insert_size") <= 380 and not self.option("its_primer"):
            # fq1_path = self.tools[0].option("out_fq").prop["path"]
            # fq2_path = self.tools[1].option("out_fq").prop["path"]
            tool_fq0 = self.tools[0].option("out_fq").prop["path"]
            tool_fq1 = self.tools[1].option("out_fq").prop["path"]
            if tool_fq0.endswith("trim.1.fq"):
                fq1_path = tool_fq0
            else:
                fq2_path = tool_fq0
            if tool_fq1.endswith("trim.2.fq"):
                fq2_path = tool_fq1
            else:
                fq1_path = tool_fq1
        else:
            for f in os.listdir(self.fastp.output_dir):
                if f.endswith("trim.1.fq"):
                    fq1_path = os.path.join(self.fastp.output_dir, f)
                if f.endswith("trim.2.fq"):
                    fq2_path = os.path.join(self.fastp.output_dir, f)
        if self.option("its_primer"):
            options["max_lenth"] = "300"
        if os.path.getsize(fq1_path) == 0:
            self.end()
        elif os.path.getsize(fq2_path) == 0:
            self.end()
        else:
            options["fq1"] = fq1_path
            options["fq2"] = fq2_path
            self.flash = self.add_tool("datasplit_v2.flash")
            self.flash.set_options(options)
            self.flash.on("end", self.run_split_by_barcode)
            self.flash.on("end", self.set_output, "flash")
            self.flash.run()

    def run_split_by_barcode(self):  # 拼接，拆样本
        options = {
            "fq": self.flash.option("out_fq"),
            "barcode_info": self.option("primer_info").prop["path"],
            "lib_name": self.option("lib_name"),
            "lib_type": self.option("lib_type")
        }
        #增加引物错配 modify by qinwen 20211012
        if self.option("mismatch"):
            options["mismatch"] = self.option("mismatch")
        self.split_by_barcode = self.add_tool("datasplit_v2.split_by_barcode")
        self.split_by_barcode.set_options(options)
        self.split_by_barcode.on("end", self.run_trim_length)
        self.split_by_barcode.on("end", self.set_output, "split_by_barcode")
        self.split_by_barcode.run()

    def run_trim_length(self):  # 拼接，长度过滤
        options = {
            "fq":  self.split_by_barcode.option("out_fq"),
            "min_len": (str(self.option("min_len")) if self.option("min_len") else str(self.min_length)),
            "lib_name": self.option("lib_name"),
        }
        if self.option("its_primer"):
            options["min_len"] = str(self.option("valid_len")) if self.option("valid_len") else "140"
        self.trim_fq_length = self.add_tool("datasplit_v2.trim_length")
        self.trim_fq_length.set_options(options)
        self.trim_fq_length.on("end", self.run_fastq_extract)
        self.trim_fq_length.run()

    def run_fastq_extract(self):  # 拼接 拆样本，统计样本的数据量
        self.logger.info(self.trim_fq_length.option("out_fq").prop["path"])
        if os.path.getsize(self.trim_fq_length.option("out_fq").prop["path"]) == 0:
            self.end()
        else:
            options = {
                "in_fastq": self.trim_fq_length.option("out_fq").prop["path"],
                "sample_primer": self.option("sample_primer"),
            }
            self.fastq_extract = self.add_tool("datasplit_v2.fastq_extract")
            self.fastq_extract.set_options(options)
            self.fastq_extract.on("end", self.set_output, "fastq_extract")
            self.fastq_extract.run()

    def run_R1_split_by_barcode(self):  # R1,拆样本
        self.logger.info("进行单端拆样本")
        options = {
            "barcode_info": self.option("primer_info").prop["path"],
            "lib_name": self.option("lib_name"),
            "split_type": "Single",
            "lib_type": self.option("lib_type")
        }
        for f in os.listdir(self.fastp.output_dir):
            if f.endswith("trim.1.fq"):
                options["fq"] = os.path.join(self.fastp.output_dir, f)
        #增加引物错配 modify by qinwen 20211012
        if self.option("mismatch"):
            options["mismatch"] = self.option("mismatch")
        self.R1_split_by_barcode = self.add_tool("datasplit_v2.split_by_barcode")
        self.R1_split_by_barcode.set_options(options)
        self.R1_split_by_barcode.on("end", self.run_R1_trim_length)
        self.R1_split_by_barcode.on("end", self.set_output, "split_by_barcode")
        self.R1_split_by_barcode.run()

    def run_R1_trim_length(self):  # 长度过滤
        options = {
            "fq":  self.R1_split_by_barcode.option("out_fq"),
            "min_len": (str(self.option("min_len")) if self.option("min_len") else str(self.min_length)),
            "lib_name": self.option("lib_name"),
        }
        self.R1_trim_length= self.add_tool("datasplit_v2.trim_length")
        self.R1_trim_length.set_options(options)
        self.R1_trim_length.on("end", self.run_R1_fastq_extract)
        self.R1_trim_length.run()

    def run_R1_fastq_extract(self):  # 拆样本，统计样本的数据量
        options = {
            "in_fastq": self.R1_trim_length.option("out_fq"),
            "sample_primer": self.option("sample_primer"),
        }
        self.R1_fastq_extract = self.add_tool("datasplit_v2.fastq_extract")
        self.R1_fastq_extract.set_options(options)
        self.R1_fastq_extract.on("end", self.set_output, "fastq_extract")
        self.R1_fastq_extract.run()

    def run(self):
        super(SingleMetaQcModule, self).run()
        self.check_primer_info()
        self.check_barcode_info()
        # if self.option("split_type") == "Auto" and self.option("lib_insert_size") > 590:
        if self.option("split_type") == "Auto" and self.option("lib_insert_size") >= 550:
            self.split_type = "Single"
        elif self.option("split_type") == "Single":
            self.split_type = "Single"
        else:
            self.split_type = "Pair"
        for fq in os.listdir(self.option("fq_dir").prop["path"]):
            if re.search(r".R1.raw.", fq):
                self.fq_list.insert(0, os.path.join(self.option("fq_dir").prop["path"], fq))
            elif fq.endswith("1.fq.gz"):
                self.fq_list.insert(0, os.path.join(self.option("fq_dir").prop["path"], fq))
            elif re.search(r".R2.raw.", fq):
                self.fq_list.append(os.path.join(self.option("fq_dir").prop["path"], fq))
            elif fq.endswith("2.fq.gz"):
                self.fq_list.append(os.path.join(self.option("fq_dir").prop["path"], fq))
            elif re.search(r"_R1_", fq):
                self.fq_list.insert(0, os.path.join(self.option("fq_dir").prop["path"], fq))
            elif re.search(r"_R2_", fq):
                self.fq_list.append(os.path.join(self.option("fq_dir").prop["path"], fq))
        if len(self.fq_list) != 2:
            raise Exception("文库文件夹中没有对应的R1,R2序列，请核实！")
        if self.option("lib_insert_size") < 200:
            self.min_length = self.option("lib_insert_size") - 20
        else:
            self.min_length = 200
        if 0 < self.option("lib_insert_size") <= 220:
            self.trim_length = 160
        elif 220 < self.option("lib_insert_size") <= 320:
            self.trim_length = 200
        elif 320 < self.option("lib_insert_size") <= 380:
            self.trim_length = 250
        if self.option("lib_type") == "official":
            self.raw_fastq_end = True
            self.run_fastp()
        else:
            self.run_getvaild_bybarcode()

    def link(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_output(self, event):
        self.logger.info("设置结果目录")
        obj = event["bind_object"]
        if event["data"] == "seq_prep":
            for f in os.listdir(obj.output_dir):
                self.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        elif event["data"] == "valid":
            for f in os.listdir(obj.output_dir):
                self.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
            for f in os.listdir(obj.work_dir):
                if f.endswith("valid.stat"):
                    self.valid_stat = os.path.join(obj.work_dir, f)
        elif event["data"] == "raw_fastq":
            raw_dir = os.path.join(self.output_dir, "meta_raw")
            if not os.path.exists(raw_dir):
                os.mkdir(raw_dir)
            for f in os.listdir(obj.output_dir):
                self.link(os.path.join(obj.output_dir, f), os.path.join(raw_dir, f))
            time.sleep(5)
            self.raw_fastq_end = True
            if self.fastq_extract_end:
                self.end()
        elif event["data"] == "fastp":
            for f in os.listdir(obj.output_dir):
                self.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
                if f.endswith("json"):
                    self.fastp_json = os.path.join(obj.output_dir, f)
        elif event["data"] == "trim":
            for f in os.listdir(obj.output_dir):
                if f.endswith("trim.hist"):
                    self.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        elif event["data"] == "flash":
            for f in os.listdir(obj.output_dir):
                self.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
            self.flash_cmd = os.path.join(obj.work_dir, "flash_cmd.o")
        elif event["data"] == "split_by_barcode":
            for f in os.listdir(obj.output_dir):
                if f.endswith("split.allLen.fq"):
                    self.split_fq = os.path.join(obj.output_dir, f)
        elif event["data"] == "fastq_extract":
            fastq_dir = os.path.join(obj.output_dir, "fastq")
            clean_dir = os.path.join(self.output_dir, "meta_clean")
            if not os.path.exists(clean_dir):
                os.mkdir(clean_dir)
            if os.path.exists(fastq_dir):  # 可能没有拆分出fastq
                for f in os.listdir(fastq_dir):
                    if f.endswith("fastq.gz"):
                        # self.link(os.path.join(fastq_dir, f), os.path.join(self.output_dir, f))
                        self.link(os.path.join(fastq_dir, f), os.path.join(clean_dir, f))
            info_path = os.path.join(obj.work_dir, "info.txt")
            new = os.path.join(self.output_dir, "info.txt")
            self.link(info_path, new)
            time.sleep(5)
            self.lib_qc_stat()
            if self.raw_fastq_end:
                self.end()
            self.fastq_extract_end = True

    def end(self):
        super(SingleMetaQcModule, self).end()

    def lib_qc_stat(self):
        self.logger.info(self.valid_stat)
        self.logger.info(self.flash_cmd)
        self.logger.info(self.fastp_json)
        self.logger.info(self.split_fq)
        raw_num, chimeric_num, raw_valid_num, chimeric_rate, raw_valid_rate = 0, 0, 0, 0, 0
        if os.path.exists(self.valid_stat):
            with open(self.valid_stat, "rb") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    raw_num = int(item[0])
                    chimeric_num = int(item[1])
                    raw_valid_num = int(item[4])
                    chimeric_rate = round(float(chimeric_num) / raw_num, 4)
                    raw_valid_rate = round(float(raw_valid_num) / raw_num, 4)
        merge_num = 0
        if os.path.exists(self.flash_cmd):
            with open(self.flash_cmd, "rb") as f:
                lines = f.readlines()
                for i in range(len(lines)-1, 0, -1):
                    if re.search("Combined reads:", lines[i]):
                        self.logger.info(lines[i])
                        num = re.findall("\d+", lines[i])
                        self.logger.info(num)
                        if len(num) > 0:
                            merge_num = int(num[0])
                        break
        trim_num, q20_rate, q30_rate = 0, 0, 0
        if os.path.exists(self.fastp_json):
            r = open(self.fastp_json, "rb")
            json_dict = json.loads(r.read())
            q_base = int(json_dict["read1_before_filtering"]["total_bases"])
            q20_1 = float(json_dict["read1_before_filtering"]["q20_bases"])
            q30_1 = float(json_dict["read1_before_filtering"]["q30_bases"])
            q20_2 = float(json_dict["read2_before_filtering"]["q20_bases"])
            q30_2 = float(json_dict["read2_before_filtering"]["q30_bases"])
            trim_num = int(json_dict["read1_after_filtering"]["total_reads"])
            q20_rate = round((q20_1 + q20_2) / (q_base * 2), 4)
            q30_rate = round((q30_1 + q30_2) / (q_base * 2), 4)
        if merge_num == 0:
            merge_num = trim_num
        split_num = 0
        if os.path.exists(self.split_fq):
            with open(self.split_fq, "rb") as f:
                lines = f.readlines()
                split_num = len(lines) / 4
        if q20_rate >= 0.95:
    	    rank = "A"
        elif q20_rate >= 0.85:
    	    rank = "B"
        elif q20_rate >= 0.75:
    	    rank = "C"
        else:
    	    rank = "D"
        trim_rate, merge_rate, split_rate, high_quality_rate = 0, 0, 0, 0
        if raw_valid_num != 0:
            trim_rate = round(float(trim_num) / raw_valid_num, 4)
            merge_rate = round(float(merge_num) / trim_num, 4)
        if merge_num != 0:
            split_rate = round(float(split_num) / merge_num, 4)
        if raw_num != 0:
            high_quality_rate = round(float(split_num) / raw_num, 4)
        qc_stat = os.path.join(self.output_dir, "lib_qc_stat.xls")
        with open(qc_stat, "wb") as w:
            w.write("#Lib\tRank\tQ20\tQ30\tRaw_pair\tchimeric\tchimeric_rate\tvalid_pair\tvalid_rate\t")
            w.write("Pair_trim\tTrim_rate\tPair_merge\tmerge_rate\tSeq_split\tSplit_rate\thighQuality_rate\n")
            w.write(self.option("lib_name") + "\t" + rank + "\t" + str(q20_rate) + "\t" + str(q30_rate) + "\t" + str(raw_num) + "\t")
            w.write(str(chimeric_num) + "\t" + str(chimeric_rate) + "\t" + str(raw_valid_num) + "\t" + str(raw_valid_rate) + "\t")
            w.write(str(trim_num) + "\t"+ str(trim_rate) + "\t" + str(merge_num) + "\t" + str(merge_rate) + "\t" + str(split_num) + "\t")
            w.write(str(split_rate) + "\t" + str(high_quality_rate) + "\n")
