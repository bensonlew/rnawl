#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
import glob
import shutil
import unittest
import re


class RnaseqMappingModule(Module):
    """
    refRNA的比对工具: tophat、hisat、star(20190506新增)
    author: shicaiping
    """
    def __init__(self, work_id):
        super(RnaseqMappingModule, self).__init__(work_id)
        options = [
            {"name": "ref_genome", "type": "string"},  # 参考基因组，在页面上呈现为下拉菜单中的选项
            {"name": "genome_version", "type": "string", "default": "Custom"},  # 参考基因组版本
            {"name": "genome_annot_version", "type": "string", "default": "Custom"},  # 参考基因组注释版本
            {"name": "mapping_method", "type": "string"},  # 比对方法，分为tophat、hisat、star
            {"name": "seq_method", "type": "string"},  # 双端测序还是单端测序
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "bam_output", "type": "outfile", "format": "align.bwa.bam_dir"},  # 输出的bam
            {"name": "bamlist", "type": "outfile", "format": "ref_rna_v2.common"},  # 输出的bamlist
            {"name": "assemble_method", "type": "string", "default": "None"},  # 拼接方法，cufflinks or stringtie or none
            {"name": "mate_std", "type": "int", "default": 50},  # 末端配对插入片段长度标准差
            {"name": "mid_dis", "type": "int", "default": 50},  # 两个成对引物间的距离中间值
            {"name": "result_reserved", "type": "int", "default": 1},  # 最多保留的比对结果数目
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {'name': 'strand_direct', 'type': 'string', 'default': 'none'},
            {"name": "strand_specific", "type": "bool", "default": False},  # 链特异性参数

            {'name': 'mapping_ratio', 'type': 'float', 'default': 60.0},
            {'name': 'mapping_sample_percent', 'type': 'float', 'default': None},
        ]
        self.add_option(options)
        self.samples = {}
        self.tool_opts = {}
        self.tools = []

    def end(self):
        n = 0.0
        for i, line in enumerate(open(os.path.join(self.output_dir, 'Comparison_results'))):
            if i and float(line.split('(')[1].split('%')[0]) < self.option('mapping_ratio'):
                n += 1
        else:
            self.option('mapping_sample_percent', n / (i + 1) * 100)
        super(RnaseqMappingModule, self).end()

    def check_options(self):
        """
        检查参数
        """
        if self.option("fastq_dir").is_set:
            self.samples = self.get_list()
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if self.option('seq_method') == "PE" and row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code = "23701201")
        if self.option('seq_method') not in ['PE', 'SE']:
            raise OptionError("请说明序列类型，PE or SE?", code = "23701203")
        if not self.option("fastq_dir").is_set and self.option('seq_method') in ["PE"]:
            if self.option("single_end_reads").is_set:
                raise OptionError("您上传的是单端测序的序列，请上传双端序列", code = "23701204")
            elif not (self.option("left_reads").is_set and self.option("right_reads").is_set):
                raise OptionError("您漏了某端序列", code = "23701205")
        if not self.option("fastq_dir").is_set and self.option('seq_method') == "SE":
            if not self.option("single_end_reads").is_set:
                raise OptionError("请上传单端序列", code = "23701206")
            elif self.option("left_reads").is_set or self.option("right_reads").is_set:
                raise OptionError("有单端的序列就够啦", code = "23701207")
        if not self.option("mapping_method").lower() in ["tophat", "hisat","star"]:
            raise OptionError("比对软件只能为tophat or hisat", code = "23701208")
        if not self.option("assemble_method").lower() in ["cufflinks", "stringtie", "None"]:
            raise OptionError("请选择拼接软件", code = "23701209")
        return True
        
    def get_opts(self):
        self.tool_opts = {
            "ref_genome": self.option("ref_genome"),
            "genome_version": self.option("genome_version"), # 参考基因组版本
            "genome_annot_version": self.option("genome_annot_version"),  # 参考基因组注释版本
            "mapping_method": self.option("mapping_method"),
            "seq_method": self.option("seq_method"),
            "assemble_method": self.option("assemble_method").lower()
        }
        if self.option("strand_specific"):
            self.tool_opts.update(
                {
                    "strand_specific": True,
                }
            )
            if self.option("mapping_method").lower() != "star":
                self.tool_opts.update(
                    {
                        'strand_direct': self.option('strand_direct'),
                    }
                )
        return True
    
    def run(self):
        super(RnaseqMappingModule, self).run()
        self.get_opts()
        if self.option("mapping_method").lower() == "tophat":
            self.logger.info("tophat开始运行")
            self.tool_run("tophat")
        elif self.option("mapping_method").lower() == "hisat":
            self.logger.info("hisat开始运行")
            self.tool_run("hisat")
        elif self.option("mapping_method").lower() == "star":
            self.logger.info("star开始运行")
            self.run_star()
        else:
            self.set_error("比对软件选择错误,程序退出", code = "23701210")
        super(RnaseqMappingModule, self).run()
        
    def tool_run(self, tool):
        if self.option("seq_method") == "PE":
            for f in self.samples:
                fq_l = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["r"])
                if tool == "star":
                    mapping_tool = self.add_tool('align.' + tool)
                else:
                    mapping_tool = self.add_tool('ref_rna_v2.' + tool)
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
                if isinstance(self.samples[f], dict):
                    fq_s = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]['s'])
                else:
                    fq_s = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f])
                mapping_tool = self.add_tool('ref_rna_v2.' + tool)
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


    def run_star(self):
        if self.option("ref_genome")=="custom_mode":
            opts={
                "ref_genome":self.option("ref_genome"),
                "ref_genome_custom":self.option("ref_genome_custom"),
                "ref_gtf_custom":self.option("ref_gtf_custom")
            }
        else:
            opts = {
            "ref_genome": self.option("ref_genome"),
            "genome_version": self.option("genome_version"),
            "genome_annot_version":self.option("genome_annot_version")

        }
        for f in self.samples:
            tool = self.add_tool("whole_transcriptome.longrna.star")
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
        self.on_rely(self.tools, self.set_output, tool)
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
        bam_list = self.output_dir + "/bamlist"
        with open(bam_list, "w") as w:
            for f in sorted(os.listdir(self.output_dir + "/bam")):
                f_path = self.output_dir + "/bam/" + f
                w.write(f_path + "\n")
        self.option("bam_output").set_path(self.output_dir + "/bam")
        self.option("bamlist").set_path(self.output_dir + "/bamlist")
        if not os.path.exists(self.output_dir + "/stat"):
            os.mkdir(self.output_dir + "/stat")
        stat_dir = self.output_dir + "/stat"
        if event["data"] == "tophat":
            for tool in self.tools:
                tophat_dir = tool.work_dir + "/tophat_out"
                stat_file = os.path.join(tophat_dir, "align_summary.txt")
                # new_file = os.path.join(stat_dir, tool.option("sample") + ".stat")
                # 更改获取sample_name的方式，避免sample_name的随机性
                sample_name = os.path.basename(glob.glob(tool.output_dir + "/*.bam")[0]).split(".bam")[0]
                new_file = os.path.join(stat_dir, sample_name + ".stat")
                if os.path.exists(new_file):
                    os.remove(new_file)
                os.link(stat_file, new_file)
        elif event["data"] == "hisat":
            self.logger.info("设置hisat输出文件")
            for tool in self.tools:
                stat_file = os.path.join(tool.work_dir, "hisat_mapping.o")
                # new_file = os.path.join(stat_dir, tool.option("sample") + ".stat")
                # 更改获取sample_name的方式，避免sample_name的随机性
                sample_name = os.path.basename(glob.glob(tool.output_dir + "/*.bam")[0]).split(".bam")[0]
                new_file = os.path.join(stat_dir, sample_name + ".stat")
                self.logger.info(stat_file)
                self.logger.info(new_file)
                if os.path.exists(new_file):
                    os.remove(new_file)
                os.link(stat_file, new_file)
        else:
             self.logger.info("设置star输出文件")
             for tool in self.tools:
                 stat_file=os.path.join(tool.work_dir, "Log.final.out")
                 # new_file = os.path.join(stat_dir, tool.option("sample") + ".stat")
                 # 更改获取sample_name的方式，避免sample_name的随机性
                 sample_name = os.path.basename(glob.glob(tool.output_dir + "/*.bam")[0]).split(".bam")[0]
                 new_file = os.path.join(stat_dir, sample_name + ".stat")
                 self.logger.info(stat_file)
                 self.logger.info(new_file)
                 if os.path.exists(new_file):
                    os.remove(new_file)
                 os.link(stat_file, new_file)

        if event["data"] == "hisat":
            self.logger.info("hisat统计结果汇总")
            Comparison_results=self.output_dir + "/Comparison_results"
            with open(Comparison_results,"w") as map_info:
             map_info.write("Sample"+"\t"+"Total_reads"+"\t"+"Total_mapped"+"\t"+"Multiple_mapped"+"\t"+"Unique_mapped"+"\n")
             stat_files=sorted(glob.glob("{}/*".format(stat_dir)))
             for fs in stat_files:
                specimen_name = os.path.basename(fs).split(".")[0]
                values = []
                f = open(fs, "r")
                for line in f:
                    if re.match(r" ", line):
                        line = line.split()
                        values.append(line[0])
                if len(values) == 13:
                    total_reads = int(values[0]) * 2
                    unmap_reads = int(values[-3])
                    map_reads = int(total_reads) - int(unmap_reads)
                    map_reads = str(map_reads) + "(" + str(float("%0.4f" % (map_reads / int(total_reads))) * 100) + "%" + ")"
                    uniq_mapped = int(values[2]) * 2 + int(values[6]) * 2 + int(values[-2])
                    uniq_mapped = str(uniq_mapped) + "(" + str(float("%0.4f" % (uniq_mapped / int(total_reads))) * 100) + "%" + ")"
                    multi_mapped = int(values[3]) * 2 + int(values[-1])
                    multi_mapped = str(multi_mapped) + "(" + str(float("%0.4f" % (multi_mapped / int(total_reads))) * 100) + "%" + ")"
                elif len(values) == 4:
                    total_reads = int(values[0])
                    unmap_reads = int(values[1])
                    map_reads = int(total_reads) - int(unmap_reads)
                    uniq_mapped = int(values[2])
                    multi_mapped = int(values[3])
                    map_reads = str(map_reads) + "(" + str(float("%0.4f" % (map_reads / int(total_reads))) * 100) + "%" + ")"
                    uniq_mapped = str(uniq_mapped) + "(" + str(float("%0.4f" % (uniq_mapped / int(total_reads))) * 100) + "%" + ")"
                    multi_mapped = str(multi_mapped) + "(" + str(float("%0.4f" % (multi_mapped / int(total_reads))) * 100) + "%" + ")"
                map_info.write(specimen_name + "\t" + str(total_reads) + "\t" + str(map_reads)+ "\t" + multi_mapped+ "\t" + uniq_mapped+ "\n")

        elif event["data"] == "tophat":
              self.logger.info("tophat统计结果汇总")
              Comparison_results = self.output_dir + "/Comparison_results"
              with open(Comparison_results, "w") as map_info:
                  map_info.write("Sample" + "\t" + "Total_reads" + "\t" + "Total_mapped" + "\t" + "Multiple_mapped" + "\t" + "Unique_mapped" + "\n")
                  stat_files = sorted(glob.glob("{}/*".format(stat_dir)))
                  for fs in stat_files:
                    specimen_name = os.path.basename(fs).split(".")[0]
                    f = open(fs, "r")
                    total_reads = 0
                    map_reads = 0
                    multiple = 0
                    for line in f:
                      if re.match(r"Left", line):
                          total_reads += int(f.next().split()[-1])
                          map_reads += int(f.next().split()[2])
                          multiple += int(f.next().split()[2])
                      if re.match(r"Right", line):
                          total_reads += int(f.next().split()[-1])
                          map_reads += int(f.next().split()[2])
                          multiple += int(f.next().split()[2])
                    total_reads_final=total_reads
                    mapping_reads_final = str(map_reads) + "(" + str(float("%0.4f" % (map_reads / total_reads)) * 100) + "%" + ")"
                    multiple_mapped_final= str(multiple) + "(" + str(float("%0.4f" % (multiple / total_reads)) * 100) + "%" + ")"
                    uniq_mapped_final = str(total_reads - multiple) + "(" + str(float("%0.4f" % ((map_reads - multiple) / total_reads)) * 100) + "%" + ")"
                    map_info.write(specimen_name +"\t"+str(total_reads)+"\t"+mapping_reads_final+"\t"+multiple_mapped_final+"\t"+uniq_mapped_final+"\n")

        else:
            self.logger.info("star统计结果汇总")
            Comparison_results = self.output_dir + "/Comparison_results"
            with open(Comparison_results, "w") as map_info:
                map_info.write("Sample"+"\t"+"Total_reads"+"\t"+"Total_mapped"+"\t"+"Multiple_mapped"+"\t"+"Unique_mapped"+"\n")
                stat_files = sorted(glob.glob("{}/*".format(stat_dir)))
                for fs in stat_files:
                    specimen_name = os.path.basename(fs).split(".")[0]
                    values = []
                    f = open(fs, "r")
                    for line in f:
                        line = line.strip()
                        if re.search("|", line):
                            line = line.split("|")
                            if len(line) > 1:
                                line_info = line[1].split()
                                line_info = line_info[0]
                                values.append(line_info)
                    if values[5] > 150:
                        total_reads = int(values[4]) * 2
                        if not values[-7].endswith("%"):
                            multi_mapped = int(values[-7]) * 2 + int(values[-9]) * 2
                        else:
                            multi_mapped = int(values[-10]) * 2 + int(values[-12]) * 2
                        uniq_mapped = int(values[6]) * 2
                        map_reads = multi_mapped + uniq_mapped
                        map_reads_rate = float("%0.4f" % (map_reads / total_reads)) * 100
                        uniq_mapped_rate = float("%0.4f" % (uniq_mapped / total_reads)) * 100
                        multi_mapped_rate = float("%0.4f" % (multi_mapped / total_reads)) * 100
                        multi_mapped_f = str(multi_mapped) + "(" + str(multi_mapped_rate) + "%" + ")"
                        map_reads_f = str(map_reads) + "(" + str(map_reads_rate) + "%" + ")"
                        uniq_mapped_f = str(uniq_mapped) + "(" + str(uniq_mapped_rate) + "%" + ")"
                        map_info.write(specimen_name + "\t" + str(total_reads) + "\t" + map_reads_f + "\t" + multi_mapped_f + "\t" + uniq_mapped_f + "\n")
                    else:
                        total_reads = int(values[4])
                        if not values[-7].endswith("%"):
                            multi_mapped = int(values[-7]) + int(values[-9])
                        else:
                            multi_mapped = int(values[-10])  + int(values[-12])
                        uniq_mapped = int(values[6])
                        map_reads = multi_mapped + uniq_mapped
                        map_reads_rate = float("%0.4f" % (map_reads / total_reads)) * 100
                        uniq_mapped_rate = float("%0.4f" % (uniq_mapped / total_reads)) * 100
                        multi_mapped_rate = float("%0.4f" % (multi_mapped / total_reads)) * 100
                        multi_mapped_f = str(multi_mapped) + "(" + str(multi_mapped_rate) + "%" + ")"
                        map_reads_f = str(map_reads)  + "(" + str(map_reads_rate) + "%" + ")"
                        uniq_mapped_f = str(uniq_mapped)  + "(" + str(uniq_mapped_rate) + "%" + ")"
                        map_info.write(specimen_name + "\t" + str(total_reads) + "\t" + map_reads_f + "\t" + multi_mapped_f + "\t" + uniq_mapped_f + "\n")
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

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "whole_transcriptome_longrna_Mapping_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "whole_transcriptome.longrna.rnaseq_mapping",
            "instant": False,
            "options": dict(
                ref_genome="Homo_sapiens",
                genome_version="GRCh38.p10",
                genome_annot_version="Ensemble_release_89",
                seq_method="PE",
                mapping_method="tophat",   # 比对软件，如star, hisat2
                fastq_dir="/mnt/ilustre/users/sanger-dev/workspace/20190926/Single_whole_transcriptome_Qc_fastp_6031/FastpRna/output/fastq",
                assemble_method="stringtie",  # 比对方法
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
        unittest.main()