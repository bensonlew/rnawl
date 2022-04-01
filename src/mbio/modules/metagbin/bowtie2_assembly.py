# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20180114
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
import shutil
from mbio.packages.metagbin.common_function import link_dir


class Bowtie2AssemblyModule(Module):
    """
    宏基因组binning组装前的比对
    """
    def __init__(self, work_id):
        super(Bowtie2AssemblyModule, self).__init__(work_id)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta_dir"},  # 输入文件,参考基因组bin_path,这是一个文件夹
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},   # 输入fastq文件夹，质控后的文件夹
            {"name": "bin_id", "type": "string"}, # bin的名称
        ]
        self.add_option(options)
        self.build_index = self.add_tool('align.bowtie_index')
        self.convert_format = self.add_module('metagbin.convert_format')
        self.choose_fastq = self.add_tool('metagbin.choose_fastq')
        self.align_tools = []
        self.cd_hit_tools = []
        self.align_result_path = ''
        self.split_fastq_tools = []
        self.cd_hit_path = ''

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option('ref_fa').is_set:
            raise OptionError('必须输入ref_fa参考序列文件夹')
        if not self.option('fastq_dir').is_set:
            raise OptionError('必须输入fastq_dir序列文件')
        if not self.option('bin_id'):
            raise OptionError('必须输入bin名称')

    def run_build_index(self):
        """
        对参考序列建index
        :return:
        """
        self.bin_id = str(self.option('bin_id')) + '.fa'
        ref_fa_path = os.path.join(self.option('ref_fa').prop['path'], self.bin_id)
        opts = ({
            'ref_fa': ref_fa_path
        })
        self.build_index.set_options(opts)
        self.build_index.on('end', self.run_split_fastq)
        self.build_index.run()

    def run_split_fastq(self):
        """
        对fastq文件进行拆分
        :return:
        """
        if os.path.exists(self.work_dir + '/all'):
            shutil.rmtree(self.work_dir + '/all')
        os.mkdir(self.work_dir + '/all')
        for file in os.listdir(self.option('fastq_dir').prop['path']):
            if re.search(r'r.fq|l.fq',file):
                fastq_path = os.path.join(self.option('fastq_dir').prop['path'], file)
                os.link(fastq_path,self.work_dir + '/all/' + file)
        for f in os.listdir(self.work_dir + '/all'):
            file_path = os.path.join(self.work_dir + '/all/', f)
            split_fastq = self.add_tool('metagbin.metagbin_split_fastq')
            opts = ({
                'in_fastq':file_path,
            })
            if f.split(".")[-2] == "r":
                opts['read_type']='2'
            elif f.split(".")[-2] == "l":
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
                new_path = os.path.join(self.work_dir, os.path.basename(file_path))
                if os.path.exists(new_path):
                    shutil.rmtree(new_path)
                link_dir(file_path, new_path)
        self.read1_path = (self.work_dir +"/fastq_1")
        self.read2_path = (self.work_dir +"/fastq_2")
        self.bin_id = str(self.option('bin_id')) + '.fa'
        ref_fasta_path = os.path.join(self.build_index.output_dir, self.bin_id)
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
        for file in os.listdir(self.option('fastq_dir').prop['path']):
            if re.search(r's.fq',file):
                file_path = os.path.join(self.option('fastq_dir').prop['path'], file)
                bowtie2 = self.add_tool('metagbin.bowtie2_align')
                opts = ({
                    'database': ref_fasta_path,
                    'fastqs': file_path,
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
        将比对结果Sam文件转为fasta和fastq文件,
        并将所有的bam文件进行合并
        :return:
        """
        sam_dir = self.align_result_path
        for i in self.align_tools:
            for f in os.listdir(i.output_dir):
                file_path = os.path.join(i.output_dir, f)
                new_path = os.path.join(self.align_result_path, os.path.basename(file_path))
                if os.path.exists(new_path):
                    os.remove(new_path)
                os.link(file_path, new_path)
        self.sample_name = self.option('bin_id') + '.fa'
        ref_path = os.path.join(self.option('ref_fa').prop['path'], self.sample_name)
        opts = ({
            "sam": sam_dir,
            "fastq_dir": self.option('fastq_dir'),
            "ref_fa": ref_path,
        })
        self.convert_format.set_options(opts)
        self.convert_format.on('end', self.run_cdhit)
        self.convert_format.run()

    def run_cdhit(self):
        """
        将fasta文件去冗余
        :return:
        """
        #sample_name = "G_" + self.option('bin_id')
        read1_path = os.path.join(self.convert_format.output_dir, "read.1.fa")
        read2_path = os.path.join(self.convert_format.output_dir, "read.2.fa")
        if os.path.basename(read1_path).split(".")[1] == "1":
            cd_hit = self.add_tool('cluster.cdhit_compare_single')
            cd_hit.set_options({
                "query": read1_path,
                "identity": 1,
                "coverage": 1,
                "pre": "fasta",
                "compare": cd_hit.work_dir,
                "qunum": 1,
            })
            self.cd_hit_tools.append(cd_hit)
        if os.path.basename(read2_path).split(".")[1] == "2":
            cd_hit = self.add_tool('cluster.cdhit_compare_single')
            cd_hit.set_options({
                "query": read2_path,
                "identity": 1,
                "coverage": 1,
                "pre": "fasta",
                "compare": cd_hit.work_dir,
                "qunum": 2,
            })
            self.cd_hit_tools.append(cd_hit)
        mk =2
        for file in os.listdir(self.option('fastq_dir').prop['path']):
            if re.search(r's.fq',file):
                mk =1
                break
        if mk == 1:
            reads_path = os.path.join(self.convert_format.output_dir, "read.s.fa")
            if os.path.basename(reads_path).split(".")[1] == "s":
                cd_hit = self.add_tool('cluster.cdhit_compare_single')
                cd_hit.set_options({
                    "query": reads_path,
                    "identity": 1,
                    "coverage": 1,
                    "pre": "fasta",
                    "compare": cd_hit.work_dir,
                    "qunum": 0,
                })
                self.cd_hit_tools.append(cd_hit)
        if len(self.cd_hit_tools) > 1:
            self.on_rely(self.cd_hit_tools, self.run_choose_fastq)
        else:
            self.cd_hit_tools[0].on('end', self.run_choose_fastq)
        for tool in self.cd_hit_tools:
            tool.run()

    def run_choose_fastq(self):
        """
        将去冗余后的fasta文件提取index，选取read1和read2相同的去提取fastq文件，
        将不同的并入single中去组装，并计算出coverage
        :return:
        """
        self.cd_hit_path = os.path.join(self.work_dir, "cd_hit")
        if os.path.exists(self.cd_hit_path):
            shutil.rmtree(self.cd_hit_path)
        os.mkdir(self.cd_hit_path)
        for i in self.cd_hit_tools:
            for f in os.listdir(i.output_dir):
                file_path = os.path.join(i.output_dir, f)
                for fr in os.listdir(file_path):
                    file_path2 = os.path.join(file_path, fr)
                    new_path = os.path.join(self.cd_hit_path, os.path.basename(file_path)+os.path.basename(file_path2))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path2, new_path)
        fasta1_path = (self.work_dir + '/cd_hit/' + 'fasta1-o')
        fasta2_path = (self.work_dir + '/cd_hit/' + 'fasta2-o')
        fastq1_path = os.path.join(self.convert_format.output_dir, "read.1.fastq")
        fastq2_path = os.path.join(self.convert_format.output_dir, "read.2.fastq")
        self.bin_id = str(self.option('bin_id')) + '.fa'
        ref_fasta_path = os.path.join(self.option('ref_fa').prop['path'], self.bin_id)
        mk =2
        for file in os.listdir(self.option('fastq_dir').prop['path']):
            if re.search(r's.fq',file):
                mk=1
                break
        if mk == 1:
            fastas_path = (self.work_dir + '/cd_hit/' + 'fasta0-o')
            fastqs_path = os.path.join(self.convert_format.output_dir, "read.s.fastq")
            self.choose_fastq.set_options({
                'fasta1': fasta1_path,
                'fasta2': fasta2_path,
                'fastas': fastas_path,
                'fastq1': fastq1_path,
                'fastq2': fastq2_path,
                'fastqs': fastqs_path,
                'ref_fa': ref_fasta_path
            })
        else:
            self.choose_fastq.set_options({
                'fasta1': fasta1_path,
                'fasta2': fasta2_path,
                'fastq1': fastq1_path,
                'fastq2': fastq2_path,
                'ref_fa': ref_fasta_path
            })
        self.choose_fastq.on('end', self.set_output)
        self.choose_fastq.run()

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info('生成结果目录')
        link_dir(self.choose_fastq.output_dir, self.output_dir)
        if os.path.exists(self.output_dir + '/Bin_coverage.xls'):
            os.remove(self.output_dir + '/Bin_coverage.xls')
        os.link(self.convert_format.output_dir + '/Bin_coverage.xls', self.output_dir + '/Bin_coverage.xls')
        self.logger.info('生成结果目录成功')
        self.end()

    def get_list(self):
        """
        根据fastq_dir路径下list.txt，将文件信息转换成字典
        dic:file_dic[sample_name][file_type]
        :return:
        """
        file_dic = dict()
        ab_rout = self.option('fastq_dir').prop['path'] + '/list.txt'
        with open(ab_rout, 'r') as list_file:
            for line in list_file:
                info = line.rstrip('\n\r').split('\t')
                name = info[1]
                type = info[2]
                if type not in ['l', 's', 'r']:
                    raise OptionError('质控样品的类型错误，必须为l/r/s之一')
                if name in file_dic.keys():
                    file_dic[name][type] = info[0]
                else:
                    file_dic[name] = {type: info[0]}
        return file_dic

    def run(self):
        super(Bowtie2AssemblyModule, self).run()
        self.run_build_index()

    def end(self):
        super(Bowtie2AssemblyModule, self).end()