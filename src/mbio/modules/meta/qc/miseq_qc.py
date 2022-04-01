# -*- coding: utf-8 -*-
# __author__ = 'xuting'

import shutil
import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class MiseqQcModule(Module):
    def __init__(self, work_id):
        super(MiseqQcModule, self).__init__(work_id)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq, sequence.fastq_dir'},  # 输入的fastq文件
            {'name': 'otu_fasta', 'type': 'outfile', 'format': 'sequence.fasta'}]  # 输出的合并到一起的fasta，供后续的otu分析用
        self.add_option(options)
        self.qc_format = self.add_tool('meta.qc.qc_format')
        self.base_info = self.add_tool('meta.qc.base_info')
        self.samples_info = self.add_tool('meta.qc.samples_info')
        self.reads_len_info = self.add_tool('meta.qc.reads_len_info')
        self.step.add_steps("fastq_format", "base_info_stat", "seq_len_stat", "sample_info_stat")

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option('in_fastq').is_set:
            raise OptionError("必须输入in_fastq参数")
        self.option('in_fastq').get_info()
        if self.get_option_object("in_fastq").format == 'sequence.fastq_dir':
            if not self.option('in_fastq').prop['has_list_file']:
                raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件')
        if self.get_option_object('in_fastq').format == 'sequence.fastq':
            if not self.option('in_fastq').prop['has_sample_info']:
                raise OptionError("fastq文件中必须在序列名中带有样本名称(以下划线分隔)")

    def qc_format_run(self):
        """
        运行qc_format，对输入的fastq文件夹或者fastq文件进行格式化
        生成fastq_dir,
        """
        self.step.fastq_format.start()
        self.step.update()
        myopt = {
            'in_fastq': self.option('in_fastq')
        }
        self.qc_format.set_options(myopt)
        self.qc_format.on("end", self.base_info_run)
        self.qc_format.on("end", self.samples_info_run)
        self.qc_format.on("end", self.reads_len_info_run)
        self.logger.info("开始对输入的fastq序列按样本进行格式化")
        self.qc_format.on('end', self.set_fasta)
        self.qc_format.run()

    def base_info_run(self):
        """
        运行碱基质量统计，对fastq文件夹里的fastq文件做碱基质量统计
        """
        self.step.fastq_format.finish()
        self.step.base_info_stat.start()
        self.step.update()
        self.logger.info("开始运行碱基质量统计")
        fastq_path = os.path.join(self.qc_format.work_dir, "output/fastq_dir")
        myopt = {
            'fastq_path': fastq_path
        }
        self.base_info.set_options(myopt)
        self.base_info.on('end', self.mv_base_info)
        self.base_info.run()

    def samples_info_run(self):
        """
        运行样品信息统计，对fasta文件夹进行样本信息统计
        """
        self.step.sample_info_stat.start()
        self.step.update()
        self.logger.info("开始进行样品信息统计")
        fasta_path = os.path.join(self.qc_format.work_dir, "output/converted_fastas")
        myopt = {
            'fasta_path': fasta_path
        }
        self.samples_info.set_options(myopt)
        self.samples_info.on('end', self.mv_samples_info)
        self.samples_info.run()

    def reads_len_info_run(self):
        """
        统计各个样本的序列长度分布
        """
        self.step.seq_len_stat.start()
        self.step.update()
        self.logger.info("开始进行序列长度分布统计")
        fasta_path = os.path.join(self.qc_format.work_dir, "output/converted_fastas")
        myopt = {
            'fasta_path': fasta_path
        }
        self.reads_len_info.set_options(myopt)
        self.reads_len_info.on('end', self.mv_reads_len)
        self.reads_len_info.run()

    def set_fasta(self):
        """
        设置合并到一起的fasta的文件路径
        """
        self.logger.info("设置fasta的文件路径")
        cat_fasta = os.path.join(self.qc_format.work_dir, "output/cat_meta.fasta")
        new_fasta = os.path.join(self.work_dir, "cat_meta.fasta")  # 因为fasta文件很大， 且涉及到尾款结算等问题，该fasta暂时不放到output中
        if os.path.exists(new_fasta):
            os.remove(new_fasta)
        os.link(cat_fasta, new_fasta)
        self.option('otu_fasta').set_path(new_fasta)

    def mv_base_info(self):
        """
        更新base_info step状态
        """
        self.step.base_info_stat.finish()
        self.step.update()

    def mv_reads_len(self):
        """
        更新reads_len_info step状态
        """
        self.step.seq_len_stat.finish()
        self.step.update()

    def mv_samples_info(self):
        """
        更新sample_info step状态
        """
        self.step.sample_info_stat.finish()
        self.step.update()

    def mv_output(self):
        output = os.path.join(self.work_dir, "output/base_info")
        if os.path.exists(output):
            shutil.rmtree(output)
        self.logger.info("开始移动baseinfo文件夹")
        base_info_dir = os.path.join(self.base_info.work_dir, "output/base_info")
        shutil.copytree(base_info_dir, output)
        output = os.path.join(self.work_dir, "output/reads_len_info")
        if os.path.exists(output):
            shutil.rmtree(output)
        self.logger.info("开始移动长度分布统计文件夹")
        reads_len_info = os.path.join(self.reads_len_info.work_dir, "output/reads_len_info")
        shutil.copytree(reads_len_info, output)
        output = os.path.join(self.work_dir, "output/samples_info")
        if os.path.exists(output):
            shutil.rmtree(output)
        self.logger.info("开始移动样品统计文件")
        samples_info_dir = os.path.join(self.samples_info.work_dir, "output/samples_info")
        shutil.copytree(samples_info_dir, output)
        self.end()

    def run(self):
        super(MiseqQcModule, self).run()
        self.qc_format_run()
        self.on_rely([self.base_info, self.samples_info, self.reads_len_info], self.mv_output)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r".", "", "结果输出目录"],
            [r'./base_info', "", "样品碱基信息统计目录"],
            ["./reads_len_info", "", "样本长度分布信息文件夹"],
            ["./samples_info", "", "样本信息结果目录"],
            ["samples_info.txt", "xls", "样本信息统计文件"]
        ])
        result_dir.add_regexp_rules([
            ["\./base_info/.+\.fastxstat\.txt$", "xls", "碱基质量统计文件"],
            ["\.reads_len_info\.txt$", "xls", "样本长度分布信息文件"]
        ])
        super(MiseqQcModule, self).end()
