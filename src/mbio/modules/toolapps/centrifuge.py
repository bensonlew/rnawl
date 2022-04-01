# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class CentrifugeModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(CentrifugeModule, self).__init__(work_id)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'qc', 'type': 'bool', 'default': False},  # 是否需要质控
            {'name': 'qc_quality', 'type': 'int', 'default': 20},  # 质控质量值标准
            {'name': 'qc_length', 'type': 'int', 'default': 30},  # 质控最短序列长度
            {'name': 'rm_host', 'type': 'bool', 'default': False},  # 是否需要去除宿主
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：E.coli ,B.taurus
            {'name': 'second_ref', 'type': 'string', 'default': ''}, #参考数据库二级下拉框数据
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'},
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
            {'name': 'clean_fq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {"name": "readlen", "type": "string", "default": "100"},
        ]
        self.add_option(options)
        self.step.add_steps('data', 'centrifuge', 'sequence', 'qc_', 'rm_host', 'clean_data')
        self.sequence = self.add_module('metagenome.reads_unzip')
        self.qc = self.add_module('meta.qc.fastp_qc')
        self.rm_host = self.add_module('meta.qc.bwa_remove_host')
        self.clean_data = self.add_module("sequence.metagbin_clean_fq")
        self.kraken = self.add_module('metagenome.centrifuge')

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('in_fastq').is_set and not self.option('clean_fq').is_set:
            raise OptionError('请输入原始序列文件夹或者质控优质序列文件夹！')
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

    def run_kraken(self):
        if self.option("qc"):
            if self.option('rm_host'):
                fastq_dir = self.rm_host.output_dir
            else:
                fastq_dir = self.qc.output_dir + '/after_qc_dir'
        else:
            fastq_dir = self.clean_data.output_dir + '/data'
        opts = {
            "fa_dir": fastq_dir,
            "readlen": int(self.option("readlen")),
        }
        self.kraken.set_options(opts)
        self.kraken.on('end', self.set_output, 'centrifuge')
        self.kraken.on("end", self.end)
        self.kraken.run()

    def run(self):
        """
        运行
        :return:
        """
        super(CentrifugeModule, self).run()
        if self.option('qc'):
            self.sequence.on('end', self.run_qc)
            if self.option('rm_host'):
                self.qc.on('end', self.run_rm_host)
                self.rm_host.on("end", self.run_kraken)
            else:
                self.qc.on('end', self.run_kraken)
            self.run_sequence()
        else:
            self.clean_data.on("end", self.run_kraken)
            self.run_clean_sequence()

    def set_output(self, event):
        """
        运行
        :return:
        """
        obj = event['bind_object']
        if event['data'] == 'centrifuge':
            if len(os.listdir(self.output_dir)) >0 :
                for i in os.listdir(self.output_dir):
                    os.remove(self.output_dir + "/" + i)
            for file in os.listdir(self.kraken.output_dir):
                os.link(self.kraken.output_dir + "/" + file, self.output_dir + "/" + file)

    def end(self):
        repaths = [
            [".", "", "Centrifuge分析结果文件目录"],
        ]
        regexps = []
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(CentrifugeModule, self).end()