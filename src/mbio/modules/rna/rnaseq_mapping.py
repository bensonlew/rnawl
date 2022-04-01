#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
import glob
import shutil
import json


class RnaseqMappingModule(Module):
    """
    refRNA的比对工具: tophat、hisat,并将用户上传的参考基因组gff文件转换为bed格式文件和gtf格式文件
    version 1.0
    author: sj
    last_modify: 2016.09.13
    """
    def __init__(self, work_id):
        super(RnaseqMappingModule, self).__init__(work_id)
        options = [
            {"name": "ref_genome", "type": "string"},  # 参考基因组，在页面上呈现为下拉菜单中的选项
            {"name": "ref_genome_custom", "type": "infile", "format": "sequence.fasta"},
            # 自定义参考基因组，用户选择customer_mode时，需要传入参考基因组
            {"name": "mapping_method", "type": "string"},  # 测序手段，分为tophat测序和hisat测序    
            {"name": "seq_method", "type": "string"},  # 双端测序还是单端测序
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "single_end_reads", "type": "infile", "format": "sequence.fastq"},  # 单端序列
            {"name": "left_reads", "type": "infile", "format": "sequence.fastq"},  # 双端测序时，左端序列
            {"name": "right_reads", "type": "infile", "format": "sequence.fastq"},  # 双端测序时，右端序列
            {"name": "bam_output", "type": "outfile", "format": "align.bwa.bam_dir"},  # 输出的bam
            {"name": "assemble_method", "type": "string", "default": "None"},  # 拼接手段，None
            {"name": "mate_std", "type": "int", "default": 50},  # 末端配对插入片段长度标准差
            {"name": "mid_dis", "type": "int", "default": 50},  # 两个成对引物间的距离中间值
            {"name": "result_reserved", "type": "int", "default": 1},  # 最多保留的比对结果数目
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "strand_specific", "type": "bool", "default": False}  # 链特异性参数
        ]
        self.add_option(options)
        self.samples = {}
        self.tool_opts = {}
        self.tools = []
        
    def check_options(self):
        """
        检查参数
        """
        if self.option("fastq_dir").is_set:
            self.samples = self.get_list()
            # self.logger.info(self.samples)
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if self.option('seq_method') == "PE" and row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列")
            elif self.option('seq_method') == "SE" and row_num != 2:
                raise OptionError("SE序列list文件应该包括文件名、样本名两列")
        if self.option('seq_method') not in ['PE', 'SE']:
            raise OptionError("请说明序列类型，PE or SE?")
        if not self.option("fastq_dir").is_set and self.option('seq_method') in ["PE"]:
            if self.option("single_end_reads").is_set:
                raise OptionError("您上传的是单端测序的序列，请上传双端序列")
            elif not (self.option("left_reads").is_set and self.option("right_reads").is_set):
                raise OptionError("您漏了某端序列")
        if not self.option("fastq_dir").is_set and self.option('seq_method') == "SE":
            if not self.option("single_end_reads").is_set:
                raise OptionError("请上传单端序列")
            elif self.option("left_reads").is_set or self.option("right_reads").is_set:
                raise OptionError("有单端的序列就够啦")
        if not self.option("mapping_method") in ["tophat", "hisat", "star"]:
            raise OptionError("tophat、hisat、star,选一个吧")
        if not self.option("assemble_method") in ["cufflinks", "stringtie", "None"]:
            raise OptionError("请选择拼接软件")
        return True
        
    def get_opts(self):
        self.tool_opts = {
            "ref_genome": self.option("ref_genome"),
            "mapping_method": self.option("mapping_method"),
            "seq_method": self.option("seq_method"),
            "assemble_method": self.option("assemble_method"),
        }
        if self.option("strand_specific"):
            self.tool_opts.update(
                {
                    "strand_specific": True
                }
            )
        return True
    
    def run(self):
        super(RnaseqMappingModule, self).run()
        self.get_opts()
        if self.option("ref_genome") == "customer_mode":
            self.tool_opts.update({
                "ref_genome_custom": self.option("ref_genome_custom")
            })
        if self.option("mapping_method") == "tophat":
            self.logger.info("tophat开始运行")
            self.tool_run("tophat")
        elif self.option("mapping_method") == "hisat":
            self.logger.info("hisat开始运行")
            self.tool_run("hisat")
        elif self.option("mapping_method") == "star":
            self.star_index = self.add_tool("align.star_index")
            self.star_index.on("end", self.run_star)
            self.run_star_index()
        else:
            raise Exception("比对软件选择错误,程序退出")
        
    def tool_run(self, tool):
        if self.option("seq_method") == "PE":
            for f in self.samples:
                fq_l = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["r"])
                mapping_tool = self.add_tool('align.' + tool)
                self.tool_opts.update({
                    'left_reads': fq_l,
                    'right_reads': fq_r,
                    'sample': f
                })
                if tool == "tophat":
                    self.tool_opts.update({
                        "mate_std": self.option("mate_std"),
                        "mid_dis": self.option("mid_dis"),
                        "result_reserved": self.option("result_reserved")
                    })
                mapping_tool.set_options(self.tool_opts)
                self.tools.append(mapping_tool)

        else:
            for f in self.samples:
                fq_s = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f])
                mapping_tool = self.add_tool('align.' + tool)
                self.tool_opts.update({
                    'single_end_reads': fq_s,
                    'sample': f
                })
                if tool == "tophat":
                    self.tool_opts.update({
                        "mate_std": self.option("mate_std"),
                        "mid_dis": self.option("mid_dis"),
                        "result_reserved": self.option("result_reserved")
                    })
                mapping_tool.set_options(self.tool_opts)
                self.tools.append(mapping_tool)
        self.on_rely(self.tools, self.set_output, tool)
        for tool in self.tools:
            tool.run()
            
    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        return samples
        
    def end(self):
        super(RnaseqMappingModule, self).end()

    def run_star(self):
        opts = {
            "ref_genome": "customer_mode",
            "ref_genome_custom": self.option("ref_genome_custom"),  # workflow调用star时统一使用customer_mode
            "ref_gtf": self.option("ref_gtf"),
            "star_index1": self.star_index.option("star_index1")
        }
        for f in self.samples:
            tool = self.add_tool("align.star")
            if self.option("seq_method") == "PE":
                fq_l = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["r"])
                opts.update({
                    "readFilesIN1": fq_l,
                    "readFilesIN2": fq_r
                })
            else:
                fq_s = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f])
                opts.update({
                    "readFilesIN": fq_s
                })
            opts.update({
                "seq_method": self.option("seq_method"),
                "sample": f
            })
            tool.set_options(opts)
            self.tools.append(tool)
        self.on_rely(self.tools, self.set_snp_output)
        for tool in self.tools:
            tool.run()

    def run_star_index(self):
        self.logger.info("开始构建star的索引")
        opts = {
            "ref_genome": self.option("ref_genome"),
            "ref_genome_custom": self.option("ref_genome_custom"),
            "ref_gtf": self.option("ref_gtf"),
            "seq_method": self.option("seq_method")
        }
        self.star_index.set_options(opts)
        self.star_index.run()
        
    def set_output(self, event):
        self.logger.info("set output")
        for f in glob.glob(r"{}/*".format(self.output_dir)):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        if not os.path.exists(self.output_dir + "/bam"):
            os.mkdir(self.output_dir + "/bam")
        new_path = self.output_dir + "/bam"
        for tool in self.tools:
            out_files = os.listdir(tool.output_dir)
            for f in out_files:
                f_path = os.path.join(tool.output_dir, f)
                target = os.path.join(new_path, f)
                if os.path.exists(target):
                    os.remove(target)
                os.link(f_path, target)
        self.option("bam_output").set_path(self.output_dir+ "/bam")
        if not os.path.exists(self.output_dir + "/stat"):
            os.mkdir(self.output_dir + "/stat")
        stat_dir = self.output_dir + "/stat"
        if event["data"] == "tophat":
            for tool in self.tools:
                tophat_dir = tool.work_dir + "/tophat_out"
                stat_file = os.path.join(tophat_dir, "align_summary.txt")
                new_file = os.path.join(stat_dir, tool.option("sample") + ".stat")
                if os.path.exists(new_file):
                    os.remove(new_file)
                os.link(stat_file, new_file)
        elif event["data"] == "hisat":
            self.logger.info("设置hisat输出文件")
            for tool in self.tools:
                stat_file = os.path.join(tool.work_dir, "hisat_mapping.o")
                new_file = os.path.join(stat_dir, tool.option("sample") + ".stat")
                self.logger.info(stat_file)
                self.logger.info(new_file)
                if os.path.exists(new_file):
                    os.remove(new_file)
                os.link(stat_file, new_file)
        self.logger.info("done")
        self.end()

    def set_snp_output(self):
        self.logger.info("set snp output")
        for f in glob.glob(r"{}/*".format(self.output_dir)):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        for tool in self.tools:
            out_files = os.listdir(tool.output_dir + "/sam")
            for f in out_files:
                f_path = os.path.join(tool.output_dir + "/sam", f)
                target = os.path.join(self.output_dir + "/sam", f)
                if not os.path.exists(self.output_dir + "/sam"):
                    os.mkdir(self.output_dir + "/sam")
                if os.path.exists(target):
                    os.remove(target)
                os.link(f_path, target)
            out_files = os.listdir(tool.output_dir + "/bam")
            for f in out_files:
                f_path = os.path.join(tool.output_dir + "/bam", f)
                target = os.path.join(self.output_dir + "/bam", f)
                if not os.path.exists(self.output_dir + "/bam"):
                    os.mkdir(self.output_dir + "/bam")
                if os.path.exists(target):
                    os.remove(target)
                os.link(f_path, target)
        self.option("bam_output").set_path(self.output_dir + "/bam")
        self.logger.info("设置snp的输出结束")
        self.end()