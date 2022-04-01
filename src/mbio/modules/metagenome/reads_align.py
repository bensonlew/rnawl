# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import get_info,link_dir


class ReadsAlignModule(Module):
    """
    将reads的fastq格式转化成fasta格式，再根据不同的方法进行比对
    """
    def __init__(self, work_id):
        super(ReadsAlignModule, self).__init__(work_id)
        options = [
            {"name": "reads_fq", "type": "infile", "format": "sequence.fastq"},#指控后的reads序列文件
            {"name": "method", "type": "string", "default": "blast"},  # blast或者diamond
            {'name': 'blast', 'type': 'string', 'default': 'blastn'},  # blastn,blastp,blastx
            {'name': 'database', 'type': 'string', 'default': 'NR'},  # NT,NR,Swiss-Prot,Custom
            {'name': 'ref_dir', 'type': 'infile', 'format': 'sequence.fasta_dir'},  # Custom时输入文件
            {'name': 'top_num', 'type': 'string', 'default': '1'},
            {'name': 'align_len', 'type': 'string', 'default': '50'},
            {'name': 'identity', 'type': 'string', 'default': '30'},
            {'name': 'evalue', 'type': 'string', 'default': '1e-5'}
        ]
        self.modules = []
        self.fqtofa = self.add_tool("metagenomic.fq_to_fa")
        self.blast = self.add_tool("metagenomic.blast")
        self.diamond = self.add_tool("metagenomic.diamond")
        self.merge_fa = self.add_tool("metagbin.cat_file")
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('reads_fq').is_set:
            raise OptionError('必须输入reads_fq文件！')
        if not self.option("blast") in ["blastn", "blastx"]:
            raise OptionError('必须输blast的比对方法不是blastn或者blastx！')
        self.logger.info(self.option("database"))

    def run_fqtofa(self):
        opts = {
            "fastq": self.option("reads_fq"),
        }
        self.fqtofa.set_options(opts)
        self.fqtofa.run()

    def run_merge_fa(self):
        opts = {
            "fa_dir": self.option('ref_dir'),
        }
        self.merge_fa.set_options(opts)
        self.merge_fa.run()

    def run_blast(self):
        opts = {
            "query": self.fqtofa.option("fasta"),
            "query_type": "nucl",
            "blast": self.option("blast"),
            "database": self.option("database"),
            "top_num": self.option("top_num"),
            "align_len": self.option("align_len"),
            "identity": self.option("identity"),
            "evalue": self.option("evalue"),
        }
        if self.option("database") in ['Custom']:
            opts['ref'] = self.merge_fa.output_dir + "/all.scaf.fa"
        if self.option("blast") in ['blastn']:
            opts['reference_type'] = "nucl"
        elif self.option("blast") in ['blastx']:
            opts['reference_type'] = "prot"
        self.blast.set_options(opts)
        self.blast.run()

    def run_diamond(self):
        opts = {
            "query": self.fqtofa.option("fasta"),
            "query_type": "nucl",
            "blast": self.option("blast"),
            "database": self.option("database"),
            "top_num": self.option("top_num"),
            "align_len": self.option("align_len"),
            "identity": self.option("identity"),
            "evalue": self.option("evalue"),
            "reference_type":"prot"
        }
        if self.option("database") in ['Custom']:
            opts['ref'] = self.merge_fa.output_dir + "/all.scaf.fa"
        self.diamond.set_options(opts)
        self.diamond.run()

    def run(self):
        """
        运行
        :return:
        """
        super(ReadsAlignModule, self).run()
        if self.option("method") == "blast":
            if self.option("database") in ['Custom']:
                self.fqtofa.on("end", self.run_merge_fa)
                self.merge_fa.on("end", self.run_blast)
                self.blast.on("end", self.set_output)
            else:
                self.fqtofa.on("end",self.run_blast)
                self.blast.on("end", self.set_output)
        elif self.option("method") == "diamond":
            if self.option("database") in ['Custom']:
                self.fqtofa.on("end", self.run_merge_fa)
                self.merge_fa.on("end", self.run_diamond)
                self.diamond.on("end", self.set_output)
            else:
                self.fqtofa.on("end", self.run_diamond)
                self.diamond.on("end", self.set_output)
        self.run_fqtofa()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(os.listdir(self.output_dir)) >= 1:
            for i in  os.listdir(self.output_dir):
                os.remove(self.output_dir + "/" +i)
        if self.option("method") == "blast":
            link_dir(self.blast.output_dir, self.output_dir)
        elif self.option("method") == "diamond":
            link_dir(self.diamond.output_dir, self.output_dir)
        self.end()

    def end(self):
        super(ReadsAlignModule, self).end()