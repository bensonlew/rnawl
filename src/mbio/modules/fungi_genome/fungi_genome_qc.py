#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile


class FungiGenomeQcModule(Module):
    """
    微生物基因组质控，二代数据：PE文库或MP文库；三代数据：pacbio数据
    last_modify: 2018.03.19
    """

    def __init__(self, work_id):
        super(FungiGenomeQcModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
            {"name": "pacbio_dir", "type": "infile", "format": "bacgenome.pacbio_dir"},#
            {"name": "analysis", "type": "string", "default": "uncomplete"},  ###流程分析模式complete，uncomplete
            {"name": "sample_name", "type": "string"},  # 样品名称
            {'name': 'qc_tool', 'type': 'string', 'default': 'fastp'},  # 质控的流程判断
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
        ]
        self.step.add_steps('raw_stat', 'pacbio_qc', 'high_sequence','clean_stat','fastp')
        self.add_option(options)
        self.raw_stat = self.add_module('bacgenome.raw_stat')  # 将二代数据解压、合并统计
        self.pacbio_qc = self.add_tool('fungi_genome.pacbio_clean')  #三代数据处理
        self.high_seq_qc = self.add_module('bacgenome.bac_high_seq_qc')  # 将二代数据处理
        self.clean_stat = self.add_module('fungi_genome.clean_stat')  # 将二代数据处理
        self.fastp = self.add_module('bacgenome.bacgenome_fastp')  # 二代数据fastp处理
        self.list = [self.raw_stat, self.clean_stat]
        self.list1 = [self.raw_stat, self.clean_stat, self.pacbio_qc]

    def check_options(self):
        """
        检查参数
        """
        if not self.option("analysis"):
            raise OptionError("请提供流程分析类型", code="22100901")
        if self.option("analysis") not in ['complete','Complete','uncomplete','Uncomplete']:
            raise OptionError("请提供正确的流程分析类型", code="22100902")


    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()


    def run_raw_stat(self):
        """
        计算二代原始数据统计
        :return:
        """
        self.logger.info("正在对二代原始数据统计开始")
        if self.option("fastq_dir").is_set:
            raw_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            self.raw_stat.set_options({
                'list_file': raw_path
            })
            self.raw_stat.on('end', self.set_output, 'raw_stat')
            self.raw_stat.run()
            self.step.raw_stat.finish()
            self.step.update()
            self.logger.info("二代原始数据统计结束")

    def run_pacbio_qc(self):
        """
        进行三代数据质控并统计、画图数据处理
        :return:
        """
        self.logger.info("正在对三代原始数据处理开始")
        if self.option("pacbio_dir").is_set:
            self.pacbio_qc.set_options({
                'input_fq': self.option("pacbio_dir").prop["path"] + '/all.pacbio.fq',
                'sample_name': self.option("sample_name")
            })
            self.pacbio_qc.on('end', self.set_output, 'pacbio_qc')
            self.pacbio_qc.run()
            self.step.pacbio_qc.finish()
            self.step.update()
            self.logger.info("三代原始数据处理结束")


    def run_high_sequence(self):
        self.logger.info("正在对二代原始数据处理开始")
        if self.option("fastq_dir").is_set:
            sample_info = os.path.join(self.option("fastq_dir").prop["path"], "sample_info")
            sample_list = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            self.high_seq_qc.set_options({
                'sample_path': sample_list,
                'sample_info': sample_info
            })
            self.high_seq_qc.on("end", self.set_output, 'high_sequence')
            self.high_seq_qc.run()
            self.step.high_sequence.finish()
            self.step.update()
            self.logger.info("正在对二代原始数据处理结束")

    def run_fastp(self):
        self.logger.info("正在对二代原始数据fastp处理开始")
        #sample_path = self.option("fastq_dir").prop["path"]
        sample_info = os.path.join(self.option("fastq_dir").prop["path"], "sample_info")
        self.fastp.set_options({
            'sample_path': self.option("fastq_dir"),
            'sample_info': sample_info
        })
        self.fastp.on("end", self.set_output, 'fastp')
        self.fastp.run()
        self.step.fastp.finish()
        self.step.clean_stat.start()
        self.step.update()
        self.logger.info("正在对二代原始数据fastp处理结束")

    def run_clean_stat(self):
        """
        计算二代原始数据统计
        :return:
        """
        self.logger.info("正在对二代clean数据统计开始")
        raw_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        if self.option("qc_tool") == "fastp":
            self.clean_stat.set_options({
                'raw_list': raw_path,
                'clean_list': self.fastp.option('clean_list'),
                'sample_name': self.option("sample_name"),
            })
        else:
            self.clean_stat.set_options({
                'raw_list': raw_path,
                'clean_list': self.high_seq_qc.option('clean_list'),
                'sample_name': self.option("sample_name"),
            })
        self.clean_stat.on("end",self.set_output,'clean_stat')
        self.clean_stat.run()
        self.step.clean_stat.finish()
        self.step.update()
        self.logger.info("正在对二代clean数据统计结束")

    def run(self):
        super(FungiGenomeQcModule, self).run()
        if self.option("analysis") in ['uncomplete', 'Uncomplete']:
            if self.option("qc_tool") == "fastp":
                if self.option("fastq_dir").is_set and not self.option("pacbio_dir").is_set:
                    self.on_rely(self.list, self.end)
                    self.fastp.on('end', self.run_clean_stat)
                    self.run_fastp()
                    self.run_raw_stat()
                elif self.option("pacbio_dir").is_set and not self.option("fastq_dir").is_set:
                    self.pacbio_qc.on('end', self.set_output)
                    self.run_pacbio_qc()
                elif self.option("fastq_dir").is_set and self.option("pacbio_dir").is_set:
                    self.on_rely(self.list1, self.end)
                    self.fastp.on('end', self.run_clean_stat)
                    self.run_fastp()
                    self.run_pacbio_qc()
                    self.run_raw_stat()
            else:
                if self.option("fastq_dir").is_set and not self.option("pacbio_dir").is_set:
                    self.on_rely(self.list, self.end)
                    self.high_seq_qc.on('end', self.run_clean_stat)
                    self.run_high_sequence()
                    self.run_raw_stat()
                elif self.option("pacbio_dir").is_set and not self.option("fastq_dir").is_set:
                    self.pacbio_qc.on('end', self.set_output)
                    self.run_pacbio_qc()
                elif self.option("fastq_dir").is_set and self.option("pacbio_dir").is_set:
                    self.on_rely(self.list1, self.end)
                    self.high_seq_qc.on('end', self.run_clean_stat)
                    self.run_high_sequence()
                    self.run_pacbio_qc()
                    self.run_raw_stat()

    def set_output(self, event):
        if self.option("analysis") in ['uncomplete', 'Uncomplete']:
            if self.option("fastq_dir").is_set and not self.option("pacbio_dir").is_set:
                self.logger.info("设置结果目录")
                if event['data'] == 'raw_stat':
                    self.linkdir(self.raw_stat.output_dir + '/fastx',"fastx")
                if event['data'] == 'high_sequence':
                    self.linkdir(self.high_seq_qc.output_dir,"cleandata")
                if event['data'] == 'clean_stat':
                    self.linkdir(self.clean_stat.output_dir+ '/data_QC', "data_QC")
                    self.linkdir(self.clean_stat.output_dir + '/fastx', "fastx")
                if event['data'] == 'fastp':
                    self.linkdir(self.fastp.output_dir, "cleandata")
            elif self.option("pacbio_dir").is_set and not self.option("fastq_dir").is_set:
                if event['data'] == 'pacbio_qc':
                    if not os.path.exists(self.output_dir + '/data_QC'):
                        os.mkdir(self.output_dir + '/data_QC')
                    if os.path.exists(self.output_dir + '/data_QC/' + self.option("sample_name") + '_PacBio_statistics.xls'):
                        os.remove(self.output_dir + '/data_QC/' + self.option("sample_name") + '_PacBio_statistics.xls')
                    os.link(self.pacbio_qc.output_dir + "/" + self.option("sample_name") + '.PacBio_statistics.xls' ,self.output_dir + '/data_QC/' + self.option("sample_name") + '_PacBio_statistics.xls' )
                    if not os.path.exists(self.output_dir + '/fastx'):
                        os.mkdir(self.output_dir + '/fastx')
                    if os.path.exists(self.output_dir + '/fastx/' + self.option("sample_name") + '.clean.len.xls'):
                        os.remove(self.output_dir + '/fastx/' + self.option("sample_name") + '.clean.len.xls')
                    os.link(self.pacbio_qc.output_dir + "/" + self.option("sample_name") + '.clean.len.xls',
                            self.output_dir + '/fastx/' + self.option("sample_name") + '.clean.len.xls')
                    self.end()
            elif self.option("fastq_dir").is_set and self.option("pacbio_dir").is_set:
                if event['data'] == 'raw_stat':
                    self.linkdir(self.raw_stat.output_dir + '/fastx',"fastx")
                if event['data'] == 'high_sequence':
                    self.linkdir(self.high_seq_qc.output_dir,"cleandata")
                if event['data'] == 'clean_stat':
                    self.linkdir(self.clean_stat.output_dir+ '/data_QC', "data_QC")
                    self.linkdir(self.clean_stat.output_dir + '/fastx', "fastx")
                if event['data'] == 'fastp':
                    self.linkdir(self.fastp.output_dir, "cleandata")
                if event['data'] == 'pacbio_qc':
                    if not os.path.exists(self.output_dir + '/data_QC'):
                        os.mkdir(self.output_dir + '/data_QC')
                    if os.path.exists(self.output_dir + '/data_QC/' + self.option("sample_name") + '_PacBio_statistics.xls'):  #guanqing.zou 20180904
                        os.remove(self.output_dir + '/data_QC/' + self.option("sample_name") + '_PacBio_statistics.xls')
                    os.link(self.pacbio_qc.output_dir + "/" + self.option("sample_name") + '.PacBio_statistics.xls' ,self.output_dir + '/data_QC/' + self.option("sample_name") + '_PacBio_statistics.xls' )
                    if not os.path.exists(self.output_dir + '/fastx'):
                        os.mkdir(self.output_dir + '/fastx')
                    if os.path.exists(self.output_dir + '/fastx/' + self.option("sample_name") + '.clean.len.xls'):   #guanqing.zou 20180904
                        os.remove(self.output_dir + '/fastx/' + self.option("sample_name") + '.clean.len.xls')
                    os.link(self.pacbio_qc.output_dir + "/" + self.option("sample_name") + '.clean.len.xls',
                            self.output_dir + '/fastx/' + self.option("sample_name") + '.clean.len.xls')

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def end(self):
        super(FungiGenomeQcModule, self).end()
