#-*- coding: utf-8 -*-

import os,re
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.file_sample import FileSampleFile
from mbio.packages.metagbin.common_function import link_dir

class FastpQcModule(Module):
    def __init__(self, work_id):
        super(FastpQcModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE OR SE
            {"name": "qualified_quality_phred", "type": "string", "default": "20"},
            # -q,一个碱基合格的质量值,默认表示phred质量> = Q是合格的。
            {"name": "length_required", "type": "string", "default": "50"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {"name": "cut_by_quality5", "type": "string", "default": "20"},  # -5,根据前面(5 ")的质量，允许每个读切割，默认是禁用的
            {"name": "cut_by_quality3", "type": "string", "default": "20"},  # -3,根据后面(3 ")的质量，允许每个读切割，默认是禁用的
            {"name": "cut_mean_quality", "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {"name": "n_base_limit", "type": "string", "default": "0"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {"name": "compression", "type": "string", "default": ""},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {"name": "thread", "type": "string", "default": "8"},  # -w,线程数
            {"name": "sickle_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # 设置结果文件后面要用
        ]
        self.add_option(options)
        #self.qc = self.add_tool("meta.fastp_info")  # 质控前后序列信息统计
        self.qc_tools = []
        self.end_times = 0

    def check_options(self):
        if not self.option("fastq_dir").is_set:
            raise OptionError("必须设置参数fastq_dir")
        return True

    def run_qc(self):
        """
        将文件夹中的文件按样本名称去质控
        :return:
        """
        self.samples = self.get_list()
        self.logger.info("{}".format(self.samples))
        fastq_path = self.option("fastq_dir").prop['path']
        for key in self.samples.keys():
            fq_l = os.path.join(fastq_path, self.samples[key]["l"])
            fq_r = os.path.join(fastq_path, self.samples[key]["r"])
            self.qc = self.add_tool("meta.qc.fastp")
            opts = ({
                "fq1": fq_l,
                "fq2": fq_r,
                "qualified_quality_phred": self.option("qualified_quality_phred"),
                "length_required": self.option("length_required"),
                "cut_mean_quality": self.option("cut_mean_quality"),
                "n_base_limit": self.option("n_base_limit"),
                "thread": self.option("thread"),
                "cut_by_quality3": self.option("cut_by_quality3"),
                "cut_by_quality5": self.option("cut_by_quality5"),
                "sample_name": key,
                "compression": self.option("compression")
            })
            self.qc.set_options(opts)
            self.qc_tools.append(self.qc)
        if len(self.qc_tools)> 1:
            self.on_rely(self.qc_tools, self.set_output)
        else:
            self.qc_tools[0].on('end', self.set_output)
        if self.qc_tools:
            for tool in self.qc_tools:
                tool.run()
        else:
            self.set_error("tool 队列是空的，无法进行后续的分析，流程终止")

    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        return samples

    def set_output(self, event):
        """
        设置结果目录
        :return:
        """
        if os.path.exists(os.path.join(self.output_dir, "qc_stat")):
            shutil.rmtree(os.path.join(self.output_dir, "qc_stat"))
        os.mkdir(os.path.join(self.output_dir, "qc_stat"))
        json_dir = os.path.join(self.output_dir, "qc_stat")
        if os.path.exists(os.path.join(self.output_dir, "after_qc_dir")):
            shutil.rmtree(os.path.join(self.output_dir, "after_qc_dir"))
        os.mkdir(os.path.join(self.output_dir, "after_qc_dir"))
        fastq_dir = os.path.join(self.output_dir, "after_qc_dir")
        with open(os.path.join(fastq_dir, "list.txt"), "w") as w:
            for obj in self.qc_tools:
                for file in os.listdir(obj.output_dir):
                    s = file.split('.clean')[0]
                    json_path = obj.work_dir + '/' + s + '.json'
                    self.link(json_path, "output/qc_stat/")
                    self.link(os.path.join(obj.output_dir, file), "output/after_qc_dir/")
                    if ".clean.1" in file:
                        di = 'l'
                    else:
                        di = 'r'
                    w.write(file + "\t" + s + "\t" + di + "\n")
        self.option('sickle_dir', self.output_dir + "/after_qc_dir")
        self.end()

    def run(self):
        super(FastpQcModule, self).run()
        self.run_qc()

    def end(self):
        super(FastpQcModule, self).end()
