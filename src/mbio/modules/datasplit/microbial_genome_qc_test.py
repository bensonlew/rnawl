# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict


class MicrobialGenomeQcTestModule(Module):
    """
    Fastp质控
    author: wangzhaoyue
    last_modify: 2017.12.13
    pe/mp两种文库，两种流程
    """

    def __init__(self, work_id):
        super(MicrobialGenomeQcTestModule, self).__init__(work_id)
        options = [
            {"name": "sample_path", "type": "infile", "format": "sequence.file_sample"},
            # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
            {'name': "sample_info", "type": "infile", "format": "sequence.barcode_info"},  # 样本信息表
            {'name': 'flag', 'type': "string", "default": "4"},  # 提取没有比对上的reads,此处固定取值4，软件默认0
            {'name': 'readl', "type": "string"},  # 切除序列的阈值
            {'name': 'illuminaclip', 'type': "string", "default": "2:30:10"},  # 2:30:10
            {'name': 'leading', 'type': "string", "default": "3"},  # 切除首端碱基质量小于0的碱基或者N
            {'name': 'tailing', 'type': "string", "default": "3"},  # 切除末端碱基质量小于20的碱基或者N
            {'name': 'sliding_window', 'type': "string", "default": "4:15"},
            # 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {'name': 'minlen', 'type': "string", "default": "36"},  # 最低reads长度
            {"name": "seqprep_quality", "type": "string", "default": '20'},
            {"name": "seqprep_length", "type": "string", "default": '25'},
            {"name": "adapter_a", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            {"name": "adapter_b", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},
            {"name": "sickle_quality", "type": "string", "default": '20'},
            {"name": "sickle_length", "type": "string", "default": '20'},
            {"name": "qual_type", "type": "string", "default": 'sanger'},
            {'name': 'qualified_quality_phred', 'type': "string", "default": "20"},  # fastp参数
            {'name': 'length_required', "type": "string", "default": "36"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {'name': "cut_by_quality5", "type": "string", "default": "20"},  # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality3", "type": "string", "default": "3"},  # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {'name': 'cut_mean_quality', "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {'name': 'n_base_limit', "type": "string", "default": "10"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {'name': 'compression', "type": "string", "default": "6"},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {'name': 'thread', "type": "string", "default": "8"},  # -w,线程数
        ]
        self.sample_info = {}
        self.pe_samples, self.mp_samples = [], []
        self.end_times = 0
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('sample_path').is_set:
            raise OptionError('必须输入样本文件夹对应的路径信息')
        row_num = len(open(self.option("sample_path").prop['path'], "r").readline().split())
        if row_num != 3:
            raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列")
        if not self.option('sample_info').is_set:
            raise OptionError('必须输入样本信息文件，包括样本，文库类型、插入片段长度三列')
        return True

    def get_info(self):
        """
        获得样本对应的路径信息，以及样本的
        :return:
        """
        with open(self.option("sample_info").prop['path'])as f:
            lines = f.readlines()
            for line in lines[1:]:
                tmp = line.strip().split('\t')
                self.sample_info[tmp[0]] = {"lib_size": tmp[2], "path": []}
                if "PE" in tmp[1]:
                    self.pe_samples.append(tmp[0])
                else:
                    self.mp_samples.append(tmp[0])
        with open(self.option('sample_path').prop['path'])as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                try:
                    if tmp[2] == 'l':
                        self.sample_info[tmp[1]]["path"].insert(0, tmp[0])
                    else:
                        self.sample_info[tmp[1]]["path"].append(tmp[0])
                except:
                    raise Exception('需要质控的序列样本{}没有相关的样本信息，请核实！'.format(tmp[1]))
        for s in self.pe_samples:
            if len(self.sample_info[s]["path"]) > 2:
                raise Exception('需要质控的序列样本{}有重名，请改样本名或分开质控！'.format(s))
            elif len(self.sample_info[s]["path"]) < 2:
                raise Exception('样本{}对应的R1,R2序列不全,请核实！'.format(s))
        for s in self.mp_samples:
            if len(self.sample_info[s]["path"]) > 2:
                raise Exception('需要质控的序列样本{}有重名，请改样本名或分开质控！'.format(s))
            elif len(self.sample_info[s]["path"]) < 2:
                raise Exception('样本{}对应的R1,R2序列不全,请核实！'.format(s))

    def run_single_microbial_genome_qc(self):
        for s in self.mp_samples:
            opts = {
                "fq1": self.sample_info[s]["path"][0],
                "fq2": self.sample_info[s]["path"][1],
                "sample_name": s,
                "insert_size": self.sample_info[s]["lib_size"],
                "flag": self.option('flag'),
                "illuminaclip": self.option('illuminaclip'),
                "leading": self.option('leading'),
                "tailing": self.option('tailing'),
                "sliding_window": self.option('sliding_window'),
                "minlen": self.option('minlen'),
                "seqprep_quality": self.option('seqprep_quality'),
                "seqprep_length": self.option('seqprep_length'),
                "adapter_a": self.option('adapter_a'),
                "adapter_b": self.option('adapter_b'),
                "sickle_quality": self.option('sickle_quality'),
                "sickle_length": self.option('sickle_length'),
                "qual_type": self.option('qual_type'),
            }
            if self.option("readl"):
                opts.update({"readl": self.option("readl")})
            self.single_micro_qc = self.add_module("datasplit.single_microbial_genome_qc")
            self.single_micro_qc.set_options(opts)
            self.single_micro_qc.on("end", self.set_output, "mp_qc")
            self.single_micro_qc.run()

    def run_phix_filter(self):
        for s in self.pe_samples:
            opts = {
                "fq1": self.sample_info[s]["path"][0],
                "fq2": self.sample_info[s]["path"][1],
                "sample_name": s,
                "insert_size": self.sample_info[s]["lib_size"],
                "flag": self.option('flag'),
            }
            if self.option("readl"):
                opts.update({"readl": self.option("readl")})
            self.phix_filter = self.add_tool("datasplit.phix_filter")
            self.phix_filter.set_options(opts)
            self.phix_filter.on("end", self.run_fastp, s)
            self.phix_filter.run()

    def run_fastp(self, event):
        obj = event["bind_object"]
        out_fq1 = obj.output_dir + '/' + event["data"] + ".filtphix.1.result.fq"
        out_fq2 = obj.output_dir + '/' + event["data"] + ".filtphix.2.result.fq"
        opts = {
            "fq1": out_fq1,
            "fq2": out_fq2,
            "qualified_quality_phred": self.option('qualified_quality_phred'),
            "length_required": self.option('length_required'),
            "cut_mean_quality": self.option('cut_mean_quality'),
            "n_base_limit": self.option('n_base_limit'),
            "compression": self.option('compression'),
            "thread": self.option('thread'),
        }
        if self.option('cut_by_quality5'):
            opts.update({'cut_by_quality5': self.option('cut_by_quality5')})
        if self.option('cut_by_quality3'):
            opts.update({'cut_by_quality3': self.option('cut_by_quality3')})
        self.fastp = self.add_tool("datasplit.fastp")
        self.fastp.set_options(opts)
        self.fastp.set_options(opts)
        self.fastp.on("end", self.set_output, event["data"])
        self.fastp.run()

    def run(self):
        """
        运行
        :return:
        """
        self.get_info()
        time.sleep(2)
        if self.mp_samples:
            self.run_single_microbial_genome_qc()
        if self.pe_samples:
            self.run_phix_filter()
        super(MicrobialGenomeQcTestModule, self).run()

    def linkdir(self, old, new, name=None):
        for f in os.listdir(old):
            f1 = os.path.join(old, f)
            if name:
                if f.endswith("clean.1.fastq.gz"):
                    f2 = new + "/" + name + ".clean.1.fq"
                else:
                    f2 = new + "/" + name + ".clean.2.fq"
            else:
                f2 = os.path.join(new, f)
            size = os.path.getsize(f1)
            if size > 0:
                if os.path.exists(f2):
                    os.remove(f2)
                os.link(f1, f2)
                if not name:  # mp的时候压缩文件
                    os.system("gzip -f {}".format(f2))
            else:
                self.logger.info("{}文件为空".format(f1))

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        obj = event["bind_object"]
        if event["data"] == "mp_qc":
            self.linkdir(os.path.join(obj.output_dir, "sickle"), self.output_dir)
        else:
            self.linkdir(obj.output_dir, self.output_dir, event["data"])
        self.end_times += 1
        if self.end_times == len(self.mp_samples) + 2 * len(self.pe_samples):
            self.end()

    def end(self):
        super(MicrobialGenomeQcTestModule, self).end()
