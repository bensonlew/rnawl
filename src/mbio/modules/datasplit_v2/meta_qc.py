# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190128

import os
import re
import json
import time
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class MetaQcModule(Module):
    """
    多样性对多个文库进行二次拆分及质控
    """
    def __init__(self, work_id):
        super(MetaQcModule, self).__init__(work_id)
        options = [
            {"name": "lib_path", "type": "infile", "format": "datasplit.path"},  # list,存放文库文件夹对应的路径信息
            {"name": "barcode_info", "type": "infile", "format": "sequence.barcode_info"},  # 文库中样本barcode及引物信息等
            {"name": "lib_specimen_id", "type": "string"},  # 文库-样本-样本的ObjectID,便于出现相同样本名称的处理
            # {"name": "lib_insert_size", "type": "string"},  # 文库插入片段长度
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE or SE
            {"name": "leading", "type": "string", "default": "0"},  # 切除首端碱基质量小于0的碱基或者N
            {"name": "tailing", "type": "string", "default": "20"},  # 切除末端碱基质量小于20的碱基或者N
            {"name": "sliding_window", "type": "string", "default": "50:20"},
            # 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {"name": "minlen", "type": "string", "default": "50"},  # 最低reads长度
            {"name": "valid_len", "type": "string"},  # -l,长度过滤阈值
            {"name": "min_lenth", "type": "string", "default": "10"},  # -m,两个reads之间所需的最小重叠长度，以提供可靠的重叠
            {"name": "max_lenth", "type": "string", "default": "100"},  # -M,两个reads之间的最大重叠长度
            {"name": "mismatch_rate", "type": "string", "default": "0.2"},  # -x,错配和重叠长度允许的最大比率
            {"name": "pred", "type": "string", "default": "33"},  # -p,FASTQ文件中的碱基的质量值，Pred33/Pred64.
            {"name": "thread", "type": "string", "default": "6"},  # -t,线程数
            {"name": "min_len", "type": "string"},  # -m,最小长度
            {"name": "split_type", "type": "string", "default": "Auto"},  # 拆分样本序列类型 Pair or Single or Auto
            # modifed by zengjing 20200102 SeqPrep、Trimmomatic的质控换成fastp
            {"name": "length_required", "type": "string", "default": "50"},  # -l,长度过滤参数，比此值短的读取将被丢弃，默认15
            {"name": "cut_right_mean_quality", "type": "string", "default": "20"},  # cut_right的平均质量要求,默认20
            {"name": "cut_right_window_size", "type": "string", "default": "50"},  # cut_right的窗口大小，默认4
            {"name": "cut_by_quality5", "type": "string", "default": "0"},  # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {"name": "cut_by_quality3", "type": "string", "default": "20"},  #  -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            # modifed by xueqinwen 20211011 增加允许引物错配数优化
            {'name': "mismatch", "type": "string"},  # -l，允许引物错配数
        ]
        self.add_option(options)
        self.clean_stat, self.raw_stat = [], []
        self.start_times, self.end_times = 0, 0
        self.its_primer = ["ITS1F_ITS2R", "ITS3F_ITS4R", "ITSF_ITSR", "ITS3F_ITS4OFR", "1737F_2043R",
                           "1761F_2043R", "ITS1F_2043R", "ITS-J1F_ITS-J2R", "ITS1F_ITS4R", "1761F_ITS2R",
                           "ITS86F_ITS4R", "ITS2-S2F_ITS4R", "ITS5F_ITS2R"]  # ITS引物，质控前加SeqPrep
        # self.single_primer = ["Arch344F_Arch915R", "bamoA1F_bamoA2R", "amoAF_amoAR", "MLfF_MLrR", "A189F_mb661R"]  # 单端引物
        self.single_primer = ["Arch344F_Arch915R", "amoAF_amoAR", "NL1F_NL4R", "ITS1F_ITS4R", "3NDF_1132rmodR"]
        self.qc_stat = {}

    def check_options(self):
        if not self.option("lib_path"):
            raise OptionError("必须输入文库文件夹对应的路径信息")
        # if not self.option("lib_insert_size"):
        #     raise OptionError("必须输入文库插入片段长度")
        if not self.option("barcode_info"):
            raise OptionError("必须输入文库的样本信息表")
        if self.option("split_type") not in ["Pair", "Single", "Auto"]:
            raise OptionError("拆分类型只能是Pair或Single或自动拆分")

    def get_info(self):
        """
        按照文库，将对应的文库序列及文库信息分开
        """
        self.project_lib, self.lib_name, self.project_num = {}, {}, {}
        self.all_lib = []
        with open(self.option("lib_path").prop["path"])as f:
            for line in f:
                tmp = line.strip().split("\t")
                self.lib_name[tmp[0]] = tmp[1]
        f = open(self.option("lib_specimen_id"), "r")
        self.lib_specimen_id = json.loads(f.read())
        with open(self.option("barcode_info").prop["path"])as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                if tmp[2] not in self.project_lib.keys():
                    self.project_num[tmp[2]] = 0
                    self.project_lib[tmp[2]] = {}
                pro_lib = tmp[2] + "-" + tmp[1] + "-" + tmp[-1]
                if pro_lib not in self.all_lib:
                    self.project_num[tmp[2]] += 1
                    self.all_lib.append(pro_lib)
                    self.project_lib[tmp[2]][str(self.project_num[tmp[2]])] = {
                        "lib": tmp[1], "specimen": [], "info": [], "its_primer": False,
                        "se_primer": False, "lib_type": tmp[-2]}
                if tmp[0] in self.project_lib[tmp[2]][str(self.project_num[tmp[2]])]["specimen"]:
                    self.project_num[tmp[2]] += 1
                    self.project_lib[tmp[2]][str(self.project_num[tmp[2]])] = {
                        "lib": tmp[1], "specimen": [], "info": [], "its_primer": False,
                        "se_primer": False, "lib_type": tmp[-2]}
                if tmp[3] in self.its_primer and not self.project_lib[tmp[2]][str(self.project_num[tmp[2]])]["its_primer"]:
                    self.project_lib[tmp[2]][str(self.project_num[tmp[2]])]["its_primer"] = True
                if tmp[3] in self.single_primer:
                    self.project_lib[tmp[2]][str(self.project_num[tmp[2]])]["se_primer"] = True
                self.project_lib[tmp[2]][str(self.project_num[tmp[2]])]["specimen"].append(tmp[0])
                self.project_lib[tmp[2]][str(self.project_num[tmp[2]])]["info"].append(tmp)
        for project_sn in self.project_lib.keys():
            for key in self.project_lib[project_sn].keys():
                barcode_info = self.work_dir + "/" + project_sn + "--" + self.project_lib[project_sn][key]["lib"] + "." + key + ".all.barcode.txt"
                primer_info = self.work_dir + "/" + project_sn + "--" + self.project_lib[project_sn][key]["lib"] + "." + key + ".all.primer.txt"
                sample_primer = self.work_dir + "/" + project_sn + "--" + self.project_lib[project_sn][key]["lib"] + "." + key + ".all.sample_primer.json"
                self.project_lib[project_sn][key]["barcode_info"] = barcode_info
                self.project_lib[project_sn][key]["primer_info"] = primer_info
                self.project_lib[project_sn][key]["sample_primer"] = sample_primer
                with open(barcode_info, "w+") as w:
                    w.write("#Sample\tBarcode-tag\tFbarcode\tRbarcode\n")
                    for tmp in self.project_lib[project_sn][key]["info"]:
                        project_order = project_sn + ":" + tmp[-1]
                        specimen_id = self.lib_specimen_id[tmp[1]][project_order][tmp[0]+"--"+tmp[3]]
                        # sample_info = project_sn + "--" + tmp[1] + "--" + specimen_id + "--" + tmp[0] + "." + tmp[3]
                        sample_info = project_sn + "--" + tmp[1] + "--" + specimen_id + "--" + tmp[0]  # seq序列id里只能是样本名称+_+编号
                        w.write(sample_info + "\t" + tmp[8] + "\t" + tmp[9] + "\t" + tmp[11] + "\n")
                with open(primer_info, "w+") as w:
                    w.write("#Sample\tF-barcode\tLinkPrimer\tR-barcode\tReversePrimer\n")
                    for tmp in self.project_lib[project_sn][key]["info"]:
                        project_order = project_sn + ":" + tmp[-1]
                        specimen_id = self.lib_specimen_id[tmp[1]][project_order][tmp[0]+"--"+tmp[3]]
                        # sample_info = project_sn + "--" + tmp[1] + "--" + specimen_id + "--" + tmp[0] + "." + tmp[3]
                        sample_info = project_sn + "--" + tmp[1] + "--" + specimen_id + "--" + tmp[0]
                        w.write(sample_info + "\t" + tmp[9] + "\t" + tmp[10] + "\t" + tmp[11] + "\t" + tmp[12] + "\n")
                with open(sample_primer, "w+") as w:
                    primer_dict = {}
                    for tmp in self.project_lib[project_sn][key]["info"]:
                        project_order = project_sn + ":" + tmp[-1]
                        specimen_id = self.lib_specimen_id[tmp[1]][project_order][tmp[0]+"--"+tmp[3]]
                        sample_info = project_sn + "--" + tmp[1] + "--" + specimen_id + "--" + tmp[0]
                        primer_dict[sample_info] = tmp[3]
                    w.write(json.dumps(primer_dict) + "\n")

    def run_single_meta_qc(self):
        module_workdir = os.path.join(self.work_dir, "module_workdir.info")
        module_file = open(module_workdir, "w")
        for project_sn in self.project_lib.keys():
            for i in self.project_lib[project_sn].keys():
                if self.option("split_type") == "Auto" and self.project_lib[project_sn][i]["se_primer"]:
                    split_type = "Single"
                else:
                    split_type = self.option("split_type")
                self.logger.info(split_type)
                options = {
                    "fq_dir": self.lib_name[self.project_lib[project_sn][i]["lib"]],
                    "barcode_info": self.project_lib[project_sn][i]["barcode_info"],
                    "primer_info": self.project_lib[project_sn][i]["primer_info"],
                    "sample_primer": self.project_lib[project_sn][i]["sample_primer"],
                    "lib_name": self.project_lib[project_sn][i]["lib"],
                    "lib_insert_size": self.project_lib[project_sn][i]["info"][0][7],
                    "fq_type": self.option("fq_type"),
                    "leading": self.option("leading"),
                    "tailing": self.option("tailing"),
                    "sliding_window": self.option("sliding_window"),
                    "minlen": self.option("minlen"),
                    "min_lenth": self.option("min_lenth"),
                    "max_lenth": self.option("max_lenth"),
                    "mismatch_rate": self.option("mismatch_rate"),
                    "pred": self.option("pred"),
                    "thread": self.option("thread"),
                    "split_type": split_type,
                    "its_primer": self.project_lib[project_sn][i]["its_primer"],
                    "length_required": self.option("length_required"),
                    "cut_right_mean_quality": self.option("cut_right_mean_quality"),
                    "cut_right_window_size": self.option("cut_right_window_size"),
                    "cut_by_quality5": self.option("cut_by_quality5"),
                    "cut_by_quality3": self.option("cut_by_quality3")
                }
                if self.option("valid_len"):
                    options["valid_len"] = self.option("valid_len")
                #增加引物错配 modify by qinwen 20211012
                if self.option("mismatch"):
                    options["mismatch"] = self.option("mismatch")
                if self.option("min_len"):
                    options["min_len"] = self.option("min_len")
                self.logger.info(self.project_lib[project_sn][i]["lib_type"])
                if re.search("双index官方多样性文库", self.project_lib[project_sn][i]["lib_type"]):
                    options["lib_type"] = "official"
                self.single_meta_qc = self.add_module("datasplit_v2.single_meta_qc")
                self.single_meta_qc.set_options(options)
                module_file.write(project_sn + "--" + self.project_lib[project_sn][str(i)]["lib"] + "\t" + str(i) + "\t" + self.single_meta_qc.work_dir + "\n")
                self.single_meta_qc.on("end", self.set_output, project_sn + "--" + self.project_lib[project_sn][i]["lib"] + "--" + str(i))
                self.single_meta_qc.on("end", self.run_celan_stat, project_sn + "--" + self.project_lib[project_sn][i]["lib"] + "--" + str(i))
                self.single_meta_qc.on("end", self.run_raw_stat, project_sn + "--" + self.project_lib[project_sn][i]["lib"] + "--" + str(i))
                self.start_times += 1
                self.single_meta_qc.run()
        module_file.close()

    def run_raw_stat(self, event):
        """
        对每个文库拆分出来的样本进行统计
        """
        obj = event["bind_object"]
        list_file = os.path.join(self.work_dir, event["data"] + ".raw.fastq.list")
        with open(list_file, "w") as w:
            # for f in os.listdir(obj.output_dir):
            raw_dir = os.path.join(obj.output_dir, "meta_raw")
            if os.path.exists(raw_dir):
                for f in os.listdir(raw_dir):
                    if f.endswith(".raw.fastq.gz"):
                        sample, r_type = "", ""
                        if f.endswith(".R1.raw.fastq.gz"):
                            sample = f.split(".R1.raw.fastq.gz")[0]
                            r_type = "l"
                            w.write(os.path.join(raw_dir, f) + "\t" + sample + "\n")
                        if f.endswith(".R2.raw.fastq.gz"):
                            sample = f.split(".R2.raw.fastq.gz")[0]
                            r_type = "r"
                        # if sample:
                        #     w.write(os.path.join(raw_dir, f) + "\t" + sample + "\t" + r_type + "\n")
        if os.path.getsize(list_file) == 0:
            self.set_output("stat_raw:" + str(event["data"]))
        else:
            options = {
                "list_file": list_file
            }
            self.fastq_stat = self.add_tool("datasplit_v2.fastq_stat")
            self.fastq_stat.set_options(options)
            self.fastq_stat.on("end", self.set_output, "stat_raw:" + event["data"])
            self.fastq_stat.run()

    def run_celan_stat(self, event):
        """
        对每个文库拆分出来的样本进行统计
        """
        obj = event["bind_object"]
        list_file = os.path.join(self.work_dir, event["data"] + ".clean.fastq.list")
        with open(list_file, "w") as w:
            # for f in os.listdir(obj.output_dir):
            clean_dir = os.path.join(obj.output_dir, "meta_clean")
            if os.path.exists(clean_dir):
                for f in os.listdir(clean_dir):
                    if f.endswith("fastq.gz"):
                        if f.endswith("clean.1.fastq.gz"):
                            continue
                        if f.endswith("clean.2.fastq.gz"):
                            continue
                        sample_info = f.split(".fastq.gz")[0]
                        w.write(os.path.join(clean_dir, f) + "\t" + sample_info + "\n")
        if os.path.getsize(list_file) == 0:
            self.set_output("stat_clean:" + str(event["data"]))
        else:
            options = {
                "list_file": list_file
            }
            self.fastq_stat = self.add_tool("datasplit_v2.fastq_stat")
            self.fastq_stat.set_options(options)
            self.fastq_stat.on("end", self.set_output, "stat_clean:" + event["data"])
            self.fastq_stat.run()

    def run(self):
        super(MetaQcModule, self).run()
        self.get_info()
        self.run_single_meta_qc()

    def link(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def link_dir(self, old_dir, new_dir):
        if not os.path.exists(old_dir):
            os.mkdir(new_dir)
        for f1 in os.listdir(old_dir):
            old = os.path.join(old_dir, f1)
            new = os.path.join(new_dir, f1)
            self.link(old, new)
        time.sleep(5)

    def set_output(self, event):
        try:
            obj = event["bind_object"]
            self.logger.info("设置结果目录:"+event["data"])
            self.logger.info(obj.output_dir)
            if event["data"].startswith("stat_clean"):
                if os.path.exists(os.path.join(obj.output_dir, "fastq_stat.xls")):
                    with open(os.path.join(obj.output_dir, "fastq_stat.xls"), "r") as f:
                        lines = f.readlines()
                        for line in lines[1:]:
                            self.clean_stat.append(line)
            elif event["data"].startswith("stat_raw"):
                if os.path.exists(os.path.join(obj.output_dir, "fastq_stat.xls")):
                    with open(os.path.join(obj.output_dir, "fastq_stat.xls"), "r") as f:
                        lines = f.readlines()
                        for line in lines[1:]:
                            self.raw_stat.append(line)
            else:
                clean_dir = os.path.join(self.output_dir, "meta_clean")
                if not os.path.exists(clean_dir):
                    os.mkdir(clean_dir)
                raw_dir = os.path.join(self.output_dir, "meta_raw")
                if not os.path.exists(raw_dir):
                    os.mkdir(raw_dir)
                trim_hist = os.path.join(clean_dir, "trim_hist")
                if not os.path.exists(trim_hist):
                    os.mkdir(trim_hist)
                for f in os.listdir(obj.output_dir):
                    if "all.raw.valid" not in f:
                        if f.endswith("fastq.gz") or f.endswith("trim.hist"):
                            old = os.path.join(obj.output_dir, f)
                            new = os.path.join(self.output_dir, f)
                            if f.endswith("trim.hist"):
                                new = os.path.join(trim_hist, event["data"]+".trim.hist")
                            self.link(old, new)
                    if f == "meta_raw":
                        self.link_dir(os.path.join(obj.output_dir, f), raw_dir)
                    elif f == "meta_clean":
                        self.link_dir(os.path.join(obj.output_dir, f), clean_dir)
                    elif f.endswith("trim.hist"):
                        old = os.path.join(obj.output_dir, f)
                        new = os.path.join(self.output_dir, f)
                        if f.endswith("trim.hist"):
                            new = os.path.join(trim_hist, event["data"]+".trim.hist")
                        self.link(old, new)
                    time.sleep(5)
                    if f == "lib_qc_stat.xls":
                        self.lib_qc_stat(os.path.join(obj.output_dir, f))
        except:
            self.logger.info("设置结果目录:"+str(event))
            self.logger.info("样本拆分结果为空")
        self.end_times += 1
        self.logger.info(self.end_times)
        self.logger.info(self.start_times)
        if self.end_times == 3*self.start_times:
            with open(os.path.join(self.output_dir, "meta_clean/fastq_stat.xls"), "w") as w:
                w.write("#Sample_ID\tTotal_Reads\tTotal_Bases\tTotal_Reads_with_Ns\tN_Reads%\tA%\tT%\tC%\tG%\tN%\tError%\tQ20%\tQ30%\tGC%\n")
                for line in self.clean_stat:
                    w.write(line)
            with open(os.path.join(self.output_dir, "meta_clean/lib_qc_stat.xls"), "w") as w:
                w.write("#Lib\tRank\tQ20\tQ30\tRaw_pair\tchimeric\tchimeric_rate\tvalid_pair\tvalid_rate\t")
                w.write("Pair_trim\tTrim_rate\tPair_merge\tmerge_rate\tSeq_split\tSplit_rate\thighQuality_rate\n")
                for lib_name in self.qc_stat.keys():
                    item = [str(i) for i in self.qc_stat[lib_name]]
                    w.write("\t".join(item) + "\n")
            with open(os.path.join(self.output_dir, "meta_raw/fastq_stat.xls"), "w") as w:
                w.write("#Sample_ID\tTotal_Reads\tTotal_Bases\tTotal_Reads_with_Ns\tN_Reads%\tA%\tT%\tC%\tG%\tN%\tError%\tQ20%\tQ30%\tGC%\n")
                for line in self.raw_stat:
                    w.write(line)
            self.end()

    def end(self):
        super(MetaQcModule, self).end()

    def lib_qc_stat(self, qc_stat):
        with open(qc_stat, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                lib_name = item[0]
                if lib_name not in self.qc_stat.keys():
                    self.qc_stat[lib_name] = [[] for i in range(16)]
                    q20_rate = float(item[2])
                    q30_rate = float(item[3])
                    raw_num = int(item[4])
                    chimeric_num = int(item[5])
                    chimeric_rate = float(item[6])
                    raw_valid_num = int(item[7])
                    raw_valid_rate = float(item[8])
                    trim_num = int(item[9])
                    trim_rate = float(item[10])
                    merge_num = int(item[11])
                    merge_rate = float(item[12])
                    split_num = int(item[13])
                    split_rate = float(item[14])
                    high_quality_rate = float(item[15])
                else:
                    q20_rate = round((float(item[2])+self.qc_stat[lib_name][2]) / 2, 4)
                    q30_rate = round((float(item[3])+self.qc_stat[lib_name][3]) / 2, 4)
                    raw_num = int(item[4])
                    chimeric_num = int(item[5]) + self.qc_stat[lib_name][5]
                    chimeric_rate = round(float(chimeric_num) / raw_num, 4)
                    raw_valid_num = int(item[7]) + self.qc_stat[lib_name][7]
                    raw_valid_rate = round(float(raw_valid_num) / raw_num, 4)
                    trim_num = int(item[9]) + self.qc_stat[lib_name][9]
                    trim_rate = round(float(trim_num) / raw_valid_num, 4)
                    merge_num = int(item[11]) + self.qc_stat[lib_name][11]
                    merge_rate = round(float(merge_num) / trim_num, 4)
                    split_num = int(item[13]) + self.qc_stat[lib_name][13]
                    split_rate = round(float(split_num) / merge_num, 4)
                    high_quality_rate = round(float(split_num) / raw_num, 4)
                if q20_rate >= 0.95:
            	    rank = "A"
                elif q20_rate >= 0.85:
            	    rank = "B"
                elif q20_rate >= 0.75:
            	    rank = "C"
                else:
            	    rank = "D"
                self.qc_stat[lib_name][0] = lib_name
                self.qc_stat[lib_name][1] = rank
                self.qc_stat[lib_name][2] = q20_rate
                self.qc_stat[lib_name][3] = q30_rate
                self.qc_stat[lib_name][4] = raw_num
                self.qc_stat[lib_name][5] = chimeric_num
                self.qc_stat[lib_name][6] = chimeric_rate
                self.qc_stat[lib_name][7] = raw_valid_num
                self.qc_stat[lib_name][8] = raw_valid_rate
                self.qc_stat[lib_name][9] = trim_num
                self.qc_stat[lib_name][10] = trim_rate
                self.qc_stat[lib_name][11] = merge_num
                self.qc_stat[lib_name][12] = merge_rate
                self.qc_stat[lib_name][13] = split_num
                self.qc_stat[lib_name][14] = split_rate
                self.qc_stat[lib_name][15] = high_quality_rate
