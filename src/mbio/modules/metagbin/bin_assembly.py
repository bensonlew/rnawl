# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import re
import time
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
from mbio.packages.metagbin.common_function import link_dir


class BinAssemblyModule(Module):
    """
    宏基因组binning基因组组装
    """
    def __init__(self, work_id):
        super(BinAssemblyModule, self).__init__(work_id)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta_dir"},  # 输入文件,参考基因组bin_path,这是一个文件夹
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},#工作流输入的文件格式
            {"name": "bam_file", "type": "infile", "format": "metagbin.file_gz"},  # metagbin.file_gz, 交互分析输入cleandata的fastq文件夹，质控后的文件
            {"name": "bin_id", "type": "string"}, # bin的名称
            {"name": "qc_stat", "type": "infile", "format": "meta_genomic.specimen_info"},  # 输入文件，样本信息表
            {"name": "insert_size", "type": "string"},
            {"name": "max_length", "type": "string"},
            {"name": "workflow", "type": "string", "default": "workflow"},  #区分工作流和交互分析,workflow和interaction
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 组装结果文件
            {"name": "fastq1", "type": "outfile", "format": "sequence.fastq"},  # 输出文件 mapping后的fastq1序列
            {"name": "fastq2", "type": "outfile", "format": "sequence.fastq"},  # 输出文件 mapping后的fastq2序列
            {"name": "fastqs", "type": "outfile", "format": "sequence.fastq"},  # 不管有没有s端输入都会有此文件；输出文件 mapping后的fastqs序列
            {"name": "coverage", "type": "outfile", "format": "sequence.profile_table"},  # 输出文件 统计binning的测序深度
        ]
        self.add_option(options)
        self.cat = self.add_tool('metagbin.cat_reads')
        self.unzip = self.add_tool('sequence.unzip')
        self.bam_fastq = self.add_tool('metagbin.bam_fastq')

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option('ref_fa').is_set:
            raise OptionError('必须输入ref_fa参考序列文件夹')
        if self.option("workflow") == "workflow":
            if not self.option("fastq_dir").is_set:
                raise OptionError('必须输入比对序列文件或文件夹')
        else:
            if not self.option('bam_file').is_set:
                raise OptionError('必须输入比对序列文件或文件夹')
        if not self.option('bin_id'):
            raise OptionError('必须输入bin名称')

    def run_metagbin_unzip(self):
        """
        将cleandata数据进行解压，并生成fastq文件夹
        :return:
        """
        opts = ({
            "file_path": self.option("bam_file").prop['path']
        })
        self.unzip.set_options(opts)
        self.unzip.on('end', self.run_bam_fastq)
        self.unzip.run()

    def run_bam_fastq(self):
        """
        将bam文件转为fastq文件
        :return:
        """
        file = os.path.basename(self.option("bam_file").prop['path'])
        if re.search(r'.tar.gz', file):
            file_name = ".".join(file.split('.')[0:-2])
        elif re.search(r'.gz', file):
            file_name = ".".join(file.split('.')[0:-1])
        else:
            file_name = ".".join(file.split('.')[0:-1])
        bam_path = os.path.join(self.unzip.output_dir, file_name)
        opts = ({
            "bam_file": bam_path
        })
        self.bam_fastq.set_options(opts)
        self.bam_fastq.on('end', self.run_bowtie2)
        self.bam_fastq.run()

    def run_bowtie2(self):
        """
        将cleandata序列fastq进行bowtie2比对，将比对结果整理后生成的fastq文件进行整理
        :return:
        """
        if self.option("workflow") != "workflow":
            self.bowtie2 = self.add_module('metagbin.bowtie2_assembly')
            opts = ({
                'ref_fa': self.option('ref_fa'),
                "fastq_dir": self.bam_fastq.output_dir,
                "bin_id": self.option('bin_id'),
            })
        else:
            self.bowtie2 = self.add_module('metagbin.bowtie_assembly')

            opts = ({
                'ref_fa': self.option('ref_fa'),
                "fastq_dir": self.option("fastq_dir"),
                "bin_id": self.option('bin_id'),
            })
        self.bowtie2.set_options(opts)
        self.bowtie2.on('end', self.run_assembly)
        self.bowtie2.run()

    def run_assembly(self):
        """
        组装，用两款软件SPAdes和SOAPdenovo
        :return:
        """
        self.sample_name = "G_" + self.option('bin_id')
        coverage_path = self.bowtie2.output_dir + "/Bin_coverage.xls"
        read1_path = self.bowtie2.output_dir + "/last.read.1.fastq"
        read2_path = self.bowtie2.output_dir + "/last.read.2.fastq"
        reads_path = self.bowtie2.output_dir + "/last.read.s.fastq"
        max_len, ave_insert_size = self.get_dic()
        with open(coverage_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                coverage = int(line[0])
                self.logger.info("coverage:{}".format(coverage))
                self.soft = self.output_dir + "/Assembly_soft.txt"
                with open(self.soft, 'wb') as w:
                    w.write("Soft\n")
                    if coverage >= 60:
                        soft = "SOAPdenovo"
                        w.write('{}\n'.format(soft))
                    else:
                        soft = "SPAdes"
                        w.write('{}\n'.format(soft))
                if coverage >= 60:
                    self.assembly = self.add_module('metagbin.bin_soapdenovo')
                    self.assembly.set_options({
                        "fastq1": read1_path,
                        "fastq2": read2_path,
                        "fastqs": reads_path,
                        "max_rd_len": str(max_len),
                        "sample_name": self.sample_name,
                        "insert_size": str(ave_insert_size),
                    })
                else:
                    self.assembly = self.add_module('metagbin.bin_spades')
                    self.assembly.set_options({
                        "fastq1": read1_path,
                        "fastq2": read2_path,
                        "fastqs": reads_path,
                        "max_rd_len": str(max_len),
                        "sample_name": self.sample_name,
                        "insert_size": str(ave_insert_size),
                    })
                self.assembly.on('end', self.set_output)
                self.assembly.run()

    def get_dic(self):
        """
        根据reads_stat和insert_size 获得每个样本的平均读长和平均插入片段长度
        :return:
        """
        insert_size = 0
        max_size = 0
        self.config = self.output_dir + "/config.txt"
        if self.option('qc_stat').is_set:
            with open(self.option('qc_stat').prop['path']) as fr1, open(self.config, 'wb') as w:
                w.write("Insert_size\tMax_length\n")
                lines = fr1.readlines()
                line_num = len(lines[1:])
                for line in lines[1:]:
                    line_list = line.strip().split('\t')
                    insert_num = int(float(line_list[2]))
                    max_num = int(float(line_list[3]))
                    insert_size += insert_num
                    max_size += max_num
                ave_insert_size = int(float(insert_size)/line_num)
                max_len = int(float(max_size)/line_num)
                w.write("{}\t{}\n".format(ave_insert_size, max_len))
        else:
            with open(self.config, 'wb') as w:
                w.write("Insert_size\tMax_length\n")
                ave_insert_size = int(float(self.option('insert_size')))
                max_len = int(float(self.option('max_length')))
                w.write('{}\t{}\n'.format(ave_insert_size, max_len))
        return max_len, ave_insert_size

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info('正在生成结果目录')
        self.sample_name = "G_" + self.option('bin_id')
        #link_dir(self.assembly.output_dir, self.output_dir)
        if os.path.exists(self.output_dir +'/Bin_coverage.xls'):
            os.remove(self.output_dir +'/Bin_coverage.xls')
        os.link(self.bowtie2.output_dir + "/Bin_coverage.xls", self.output_dir +'/Bin_coverage.xls')
        self.option('scaffold', self.assembly.option('scaffold'))
        if os.path.exists(self.output_dir +'/'+ self.sample_name +'_scaffold.fa'):
            os.remove(self.output_dir +'/'+ self.sample_name +'_scaffold.fa')
        os.link(self.option('scaffold').prop['path'], self.output_dir +'/'+ self.sample_name +'_scaffold.fa')
        self.option('fastq1', self.bowtie2.output_dir + "/last.read.1.fastq")
        self.option('fastq2', self.bowtie2.output_dir + "/last.read.2.fastq")
        self.option('fastqs', self.bowtie2.output_dir + "/last.read.s.fastq")
        self.option('coverage', self.output_dir + "/Bin_coverage.xls")
        self.logger.info("结果目录生成成功")
        self.end()

    def run(self):
        super(BinAssemblyModule, self).run()
        if self.option("workflow") != "workflow":
            self.run_metagbin_unzip()
        else:
            self.run_bowtie2()


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BinAssemblyModule, self).end()
