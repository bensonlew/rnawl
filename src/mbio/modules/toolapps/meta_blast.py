# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.toolapps.common import link_dir

class MetaBlastModule(Module):
    """
    balst的小工具
    """
    def __init__(self, work_id):
        super(MetaBlastModule, self).__init__(work_id)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'clean_fq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'qc', 'type': 'string', 'default': "true"},  # true,false,fasta
            {'name': 'rm_host', 'type': 'bool', 'default': False},  # 是否需要去除宿主
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：E.coli ,B.taurus
            {'name': 'second_ref', 'type': 'string', 'default': ''}, #参考数据库二级下拉框数据
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'},
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
            {'name': 'fasta_type', 'type': 'string', 'default': 'nucleotide'},  # nucleotide,protein
            {'name': 'fasta_dir', 'type': 'infile', 'format': 'toolapps.fasta_dir'},
            {'name': 'blast', 'type': 'string', 'default': 'blastn'}, #blastn,blastp,blastx
            {'name': 'database', 'type': 'string', 'default': 'NT'},#NT,NR,Swiss-Prot,Custom
            {'name': 'ref_dir', 'type': 'infile', 'format': 'sequence.fasta_dir'}, #Custom时输入文件
            {'name': 'top_num', 'type': 'string', 'default': '1'},
            {'name': 'align_len', 'type': 'string', 'default': '50'},
            {'name': 'identity', 'type': 'string', 'default': '30'},
            {'name': 'evalue', 'type': 'string', 'default': '1e-5'}
        ]
        self.add_option(options)
        self.step.add_steps('data', 'blast', 'sequence', 'qc_', 'rm_host', 'clean_data', "fa_blast")
        self.sequence = self.add_module('metagenome.reads_unzip')
        self.qc = self.add_module('meta.qc.fastp_qc')
        self.rm_host = self.add_module('meta.qc.bwa_remove_host')
        self.clean_data = self.add_module("sequence.metagbin_clean_fq")

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('in_fastq').is_set and not self.option('clean_fq').is_set and not self.option('fasta_dir').is_set:
            raise OptionError('请输入原始序列文件夹或者质控优质序列文件夹或者上传fasta格式的文件夹！')
        if self.option('rm_host'):
            if self.option('ref_database') == '' and not self.option('ref_undefined').is_set:
                raise OptionError('已选择去宿主，需输入参考数据库或参考序列')
            if self.option('ref_database') not in ['', 'Custom'] and self.option('ref_undefined').is_set:
                raise OptionError('去宿主不可同时提供参考数据库及参考序列')

    def run_sequence(self):
        opts = {
            'fastq_dir': self.option('in_fastq'),
        }
        self.sequence.set_options(opts)
        self.sequence.on('end', self.set_output, 'sequence')
        self.sequence.run()

    def run_qc(self):
        opts = {
            'fastq_dir': self.sequence.output_dir + '/data',
        }
        self.qc.set_options(opts)
        self.qc.on('end', self.set_output, 'qc_')
        self.qc.run()

    def run_rm_host(self):
        ref_database = self.option('ref_database') + "," + self.option('second_ref')
        opts = {
            'fq_type': 'PE',
            'ref_database': ref_database,
            'ref_undefined': self.option('ref_undefined'),
        }
        if self.option('ref_database') == 'Custom':
            opts['ref_database'] = ""
        if self.option('qc'):
            opts['fastq_dir'] = self.qc.output_dir + "/after_qc_dir"
        else:
            opts['fastq_dir'] = self.option('in_fastq')
        self.rm_host.set_options(opts)
        self.rm_host.on('end', self.set_output, 'rm_host')
        self.rm_host.run()

    def run_clean_sequence(self):
        opts = {
            'fastq_dir': self.option('clean_fq'),
        }
        self.clean_data.set_options(opts)
        self.clean_data.on('end', self.set_output, 'clean_data')
        self.clean_data.run()

    def run_blast(self):
        fastq_dir = ""
        if self.option('qc') in ['true', 'True']:
            if self.option('rm_host'):
                fastq_dir = self.rm_host.output_dir
            else:
                fastq_dir = self.qc.output_dir + '/after_qc_dir'
        elif self.option('qc') in ['false', 'False']:
            fastq_dir = self.clean_data.output_dir + '/data'
        self.blast = self.add_module('metagenome.reads_align_database')
        opts = {
            "fa_dir": fastq_dir,
            "method": "blast",
            "blast": self.option("blast"),
            "database": self.option("database"),
            "ref_dir": self.option("ref_dir"),
            "top_num": self.option("top_num"),
            "align_len": self.option("align_len"),
            "identity": self.option("identity"),
            "evalue": self.option("evalue"),
        }
        self.blast.set_options(opts)
        self.blast.on('end', self.set_output, 'blast')
        self.blast.on("end", self.end)
        self.blast.run()

    def run_fa_blast(self):
        self.fasta_blast = self.add_module('metagenome.fa_align')
        opts = {
            "fa_dir": self.option('fasta_dir'),
            "method": "blast",
            "fasta_type": self.option("fasta_type"),
            "blast": self.option("blast"),
            "database": self.option("database"),
            "ref_dir": self.option("ref_dir"),
            "top_num": self.option("top_num"),
            "align_len": self.option("align_len"),
            "identity": self.option("identity"),
            "evalue": self.option("evalue"),
        }
        self.fasta_blast.set_options(opts)
        self.fasta_blast.on('end', self.set_output, 'fa_blast')
        self.fasta_blast.run()

    def run(self):
        """
        运行
        :return:
        """
        super(MetaBlastModule, self).run()
        if self.option('qc') in ['true', 'True']:
            self.sequence.on('end', self.run_qc)
            if self.option('rm_host'):
                self.qc.on('end', self.run_rm_host)
                self.rm_host.on("end", self.run_blast)
            else:
                self.qc.on('end', self.run_blast)
            self.run_sequence()
        elif self.option('qc') in ['false', 'False']:
            self.clean_data.on("end", self.run_blast)
            self.run_clean_sequence()
        elif self.option('qc') in ['fasta', 'Fasta']:
            self.run_fa_blast()

    def set_output(self, event):
        """
        运行
        :return:
        """
        obj = event['bind_object']
        if event['data'] == 'blast':
            if len(os.listdir(self.output_dir)) >0 :
                for i in os.listdir(self.output_dir):
                    os.remove(self.output_dir + "/" + i)
            link_dir(self.blast.output_dir, self.output_dir)
            self.end()
        if event['data'] == 'fa_blast':
            if len(os.listdir(self.output_dir)) >0 :
                for i in os.listdir(self.output_dir):
                    os.remove(self.output_dir + "/" + i)
            link_dir(self.fasta_blast.output_dir, self.output_dir)
            self.end()

    def end(self):
        repaths = [
            [".", "", "blast分析结果文件目录"],
        ]
        regexps = []
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(MetaBlastModule, self).end()