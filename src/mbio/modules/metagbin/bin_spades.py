# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'

import os
import re
import time,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
from mbio.packages.metagbin.common_function import link_dir


class BinSpadesModule(Module):
    def __init__(self, work_id):
        super(BinSpadesModule, self).__init__(work_id)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.l.fastq
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.r.fastq
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.s.fastq
            {'name': 'sample_name', "type": "string"},  # 基因组名称对应项目的genome名称
            {"name": "max_rd_len", "type": "string"},  # read最大读长
            {"name": "insert_size", "type": "string"},  # 平均插入片段长度
            {"name": "reverse_seq", "type": "string", "default": "0"},  # 配置文件的其他参数
            {"name": "asm_flags", "type": "string", "default": "3"},  # 配置文件的其他参数
            {"name": "rank", "type": "string", "default": "1"},  # 配置文件的其他参数
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
        ]
        self.add_option(options)
        self.spades_assemble = self.add_tool('metagbin.metagbin_spades')
        self.gapcloser = self.add_tool('assemble.gapcloser_scaf_bin')
        self.scaf_agp_contig = self.add_tool('assemble.scaf_agp_contig')

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fastq1').is_set:
            raise OptionError('必须输入fastq1序列')
        if not self.option('fastq2'):
            raise OptionError('必须输入fastq2序列')
        #if not self.option('fastqs'):
            #raise OptionError('必须输入fastqs序列', code="")

    def run_config(self):
        """
        配置config文件
        :return:
        """
        config_file = self.work_dir + "/" + self.option('sample_name')+ ".config"
        with open(config_file, "w+") as fw:
            first = "max_rd_len=" + self.option('max_rd_len') + "\n"
            second = "[LIB]" + "\n"
            third = "avg_ins=" + self.option('insert_size') + "\n"
            forth = "reverse_seq=" + self.option('reverse_seq') + "\n"
            fifth = "asm_flags=" + self.option('asm_flags') + "\n"
            sixth = "rank=" + self.option('rank') + "\n"
            q1 = "q1=" + self.option('fastq1').prop['path'] + "\n"
            q2 = "q2=" + self.option('fastq2').prop['path']
            if not self.option('fastqs').is_set:
                fw.write(first + second + third + forth + fifth + sixth + q1 + q2)
            else:
                qs = "\nq=" + self.option('fastqs').prop['path']
                fw.write(first + second + third + forth + fifth + sixth + q1 + q2 + qs)

    def run_assembly(self):
        """
        用SPAdes软件对序列进行组装
        :return:
        """
        if self.option('fastqs').is_set:
            opts = {
                "fastq1": self.option('fastq1'),
                "fastq2": self.option('fastq2'),
                "fastqs": self.option('fastqs'),
                "sample_name": self.option('sample_name'),
            }
        else:
            opts = {
                "fastq1": self.option('fastq1'),
                "fastq2": self.option('fastq2'),
                "sample_name": self.option('sample_name'),
            }
        self.spades_assemble.set_options(opts)
        self.spades_assemble.on('end', self.run_gapcloser)
        self.spades_assemble.run()

    def run_gapcloser(self):
        """
        对组装结果根据参考序列进行补洞
        :return:
        """
        config_file = self.work_dir + "/" + self.option('sample_name')+ ".config"
        opts = {
            "seq_scaf": self.spades_assemble.option('scaffold'),
            "config": config_file,
            "sample_name": self.option('sample_name'),
        }
        self.gapcloser.set_options(opts)
        self.gapcloser.on('end', self.run_scaf_agp_contig)
        self.gapcloser.run()

    def run_scaf_agp_contig(self):
        """
        主要是过滤掉长度小于200bp的序列，挑选最终的scaffold和contigs
        :return:
        """
        opts = {
            "seq_scaf": self.gapcloser.option('seq_scaf'),
            "sample_name": self.option('sample_name'),
        }
        self.scaf_agp_contig.set_options(opts)
        self.scaf_agp_contig.on('end', self.set_output)
        self.scaf_agp_contig.run()

    def run(self):
        """
        运行
        :return:
        """
        super(BinSpadesModule, self).run()
        self.run_config()
        self.run_assembly()

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        #if os.path.exists(self.output_dir + '/scf'):
            #shutil.rmtree(self.output_dir + '/scf')
        #link_dir(self.gapcloser.output_dir,self.output_dir + '/scf')
        #self.option('gap_scaffold').set_path(self.output_dir + '/scf/' + self.option('sample_name') + '.scaffold.fna')
        if os.path.exists(self.output_dir + '/'+self.option('sample_name')+'_scaffold.fa'):
            os.remove(self.output_dir + '/'+self.option('sample_name')+'_scaffold.fa')
        os.link(self.scaf_agp_contig.option('scaffold').prop['path'], self.output_dir + '/'+self.option('sample_name')+'_scaffold.fa')
        self.option('scaffold',self.output_dir + '/'+self.option('sample_name')+'_scaffold.fa')
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BinSpadesModule, self).end()