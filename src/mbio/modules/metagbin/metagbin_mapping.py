#-*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re,shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class MetagbinMappingModule(Module):
    """
    宏基因组binning mapping计算丰度
    """
    def __init__(self, work_id):
        super(MetagbinMappingModule, self).__init__(work_id)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 输入文件,参考基因组bin
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 输入fastq文件夹，质控后的文件夹
            {"name": "bam_dir", "type": "outfile", "format": "metagbin.bam_dir"},  #
        ]
        self.add_option(options)
        self.split_fastq_tools = []
        self.build_index = self.add_tool('align.bowtie_index')
        self.convert_format = self.add_module('metagbin.convert_format')
        self.align_tools = []
        self.convert_format_tools = []

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option('ref_fa').is_set:
            raise OptionError('必须输入ref_fa序列文件')
        if not self.option('fastq_dir').is_set:
            raise OptionError('必须输入fastq_dir序列文件')

    def run_build_index(self):
        """
        对参考序列建index
        :return:
        """
        ref_fa_path = self.option('ref_fa').prop['path']
        ref_pre = os.path.basename(self.option('ref_fa').prop['path'])
        if not os.path.exists(self.work_dir + "/ref"):
            os.mkdir(self.work_dir + "/ref/")
        if os.path.exists(self.work_dir + "/ref/" + ref_pre):
            os.remove(self.work_dir + "/ref/" + ref_pre)
        os.link(ref_fa_path, self.work_dir + "/ref/" + ref_pre)
        opts = ({
            'ref_fa': ref_fa_path,
            'bowtie_version': 'bowtie2',
        })
        self.build_index.set_options(opts)
        self.build_index.run()

    def run_split_fastq(self):
        """
        对fastq文件进行拆分
        :return:
        """
        self.samples = self.get_dic(self.option('fastq_dir').prop['path'] + '/list.txt')
        for sample, value in sorted(self.samples.items()):
            for type, path in value.items():
                split_fastq = self.add_tool('metagbin.split_fastq')
                opts = ({
                    'in_fastq': self.option('fastq_dir').prop['path'] + '/' + path,
                    'sample': sample,
                })
                if type in ['r']:
                    opts['read_type'] = '2'
                    split_fastq.set_options(opts)
                    self.split_fastq_tools.append(split_fastq)
                elif type in ['l']:
                    opts['read_type'] = '1'
                    split_fastq.set_options(opts)
                    self.split_fastq_tools.append(split_fastq)
        if len(self.split_fastq_tools) > 1:
            self.on_rely(self.split_fastq_tools, self.run_align)
        else:
            self.split_fastq_tools[0].on('end', self.run_align)
        for tool in self.split_fastq_tools:
            tool.run()

    def run_align(self):
        """
        用bowtie2软件进行比对
        :return:
        """
        self.align_result_path = os.path.join(self.work_dir, "align_result")
        if os.path.exists(self.align_result_path):
            shutil.rmtree(self.align_result_path)
        os.mkdir(self.align_result_path)
        for i in self.split_fastq_tools:
            for f in os.listdir(i.output_dir):
                file_path = os.path.join(i.output_dir, f)
                new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                link_dir(file_path, new_path)
        fa = os.path.basename(self.option('ref_fa').prop['path'])
        ref_fasta_path = self.build_index.output_dir + "/" + fa
        for sample in sorted(self.samples.keys()):
            self.read1_path = (self.align_result_path + "/" + sample + "_fastq_1")
            self.read2_path = (self.align_result_path + "/" + sample + "_fastq_2")
            for file in os.listdir(self.read1_path):
                read1_path = os.path.join(self.read1_path, file)
                read2_path = os.path.join(self.read2_path, file)
                bowtie2 = self.add_tool('metagbin.bowtie2_align')
                opts = ({
                    'database': ref_fasta_path,
                    'fastq1': read1_path,
                    'fastq2': read2_path,
                })
                bowtie2.set_options(opts)
                self.align_tools.append(bowtie2)
        for sample,value in sorted(self.samples.items()):
            for type,path in value.items():
                if type in ['s']:
                    bowtie2 = self.add_tool('metagbin.bowtie2_align')
                    opts = ({
                        'database': ref_fasta_path,
                        'fastqs': self.option('fastq_dir').prop['path'] + '/' + path,
                    })
                    bowtie2.set_options(opts)
                    self.align_tools.append(bowtie2)
        if len(self.align_tools) > 1:
            self.on_rely(self.align_tools, self.run_convert_format)
        else:
            self.align_tools[0].on('end', self.run_convert_format)
        for tool in self.align_tools:
            tool.run()

    def run_convert_format(self):
        """
        将比对结果Sam文件转为bam文件
        :return:
        """
        self.sam_path = os.path.join(self.work_dir, "sam_dir")
        if os.path.exists(self.sam_path):
            shutil.rmtree(self.sam_path)
        os.mkdir(self.sam_path)
        for i in self.align_tools:
            for f in os.listdir(i.output_dir):
                file_path = os.path.join(i.output_dir, f)
                new_path = os.path.join(self.sam_path, os.path.basename(file_path))
                if os.path.exists(new_path):
                    os.remove(new_path)
                os.link(file_path, new_path)
        opts = ({
            "sam": self.sam_path,
            "analysis": "bam",
        })
        self.convert_format.set_options(opts)
        self.convert_format.on("end", self.set_output)
        self.convert_format.run()

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info('生成结果目录')
        self.option("bam_dir", self.convert_format.output_dir + '/bam_sort')
        self.logger.info('生成结果目录成功')
        self.end()

    def run(self):
        super(MetagbinMappingModule, self).run()
        self.build_index.on("end", self.run_split_fastq)
        self.run_build_index()

    def end(self):
        super(MetagbinMappingModule, self).end()

    def get_dic(self,list_file):
        dict = {}
        with open (list_file,'r') as f:
            lines = f.readlines()
            for line in lines:
                lin = line.rstrip('\r\n').split('\t')
                if lin[1] in dict.keys():
                     if lin[2] not in dict.values():
                         dict[lin[1]][lin[2]] = lin[0]
                else:
                    dict[lin[1]] = {lin[2]: lin[0]}
        return dict


