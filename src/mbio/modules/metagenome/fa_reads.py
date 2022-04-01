# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import link_dir,add_merge


class FaReadsModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(FaReadsModule, self).__init__(work_id)
        options = [
            {"name": "in_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample", "type": "string"},
            {"name": "method", "type": "string", "default": "blast"},  # blast或者diamond
            {'name': 'fasta_type', 'type': 'string', 'default': 'nucl'},  # nucl,prot
            {'name': 'blast', 'type': 'string', 'default': 'blastn'},  # blastn,blastp,blastx
            {'name': 'database', 'type': 'string', 'default': 'NT'},  # NT,NR,Swiss-Prot,Custom
            {'name': 'ref_dir', 'type': 'infile', 'format': 'sequence.fasta_dir'},  # Custom时输入文件
            {'name': 'top_num', 'type': 'string', 'default': '1'},
            {'name': 'align_len', 'type': 'string', 'default': '50'},
            {'name': 'identity', 'type': 'string', 'default': '30'},
            {'name': 'evalue', 'type': 'string', 'default': '1e-5'},
        ]
        self.modules = []
        self.add_option(options)
        self.split_seq = self.add_tool("metagenomic.split_seq")
        self.sample_path = defaultdict(list)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('in_fasta').is_set:
            raise OptionError('必须输入in_fasta文件夹')
        if self.option('fasta_type') not in ['nucl', 'pro']:
            raise OptionError('必须输入序列类型nucleotide或protein文件！')
        else:
            if self.option('fasta_type') in ['nucl']:
                self.query_type = "nucl"
            elif self.option('fasta_type') in ['prot']:
                self.query_type = "prot"
        if self.option('blast') not in ['blastn', 'blastp', 'blastx']:
            raise OptionError('必须比对方法不是blastn，blastp, blastx其中一个！')
        else:
            if self.option('blast') in ['blastn']:
                self.reference_type = "nucl"
            elif self.option('blast') in ['blastp','blastx']:
                self.reference_type = "prot"

    def run_splitfq(self):
        opts = {
            "in_fasta": self.option('in_fasta'),
            "sample": self.option('sample'),
        }
        self.split_seq.set_options(opts)
        self.split_seq.run()
    
    def run_merge_fa(self):
        self.merge_fa = self.add_tool("metagbin.cat_file")
        opts = {
            "fa_dir": self.option('ref_dir'),
        }
        self.merge_fa.set_options(opts)
        self.merge_fa.run()
        
    def run_align(self):
        if self.num <= 500000:
            if self.option("method") in ["blast"]:
                self.blast = self.add_tool("metagenomic.blast")
                opts = {
                    "query": self.option("in_fasta"),
                    "query_type": self.query_type,
                    "blast": self.option("blast"),
                    "database": self.option("database"),
                    "reference_type": self.reference_type,
                    "top_num": self.option("top_num"),
                    "align_len": self.option("align_len"),
                    "identity": self.option("identity"),
                    "evalue": self.option("evalue"),
                }
                if self.option("database") in ['Custom']:
                    opts['ref'] = self.merge_fa.output_dir + "/all.scaf.fa"
                self.blast.set_options(opts)
                self.blast.on("end", self.set_output)
                self.blast.run()
            elif self.option("method") in ["diamond"]:
                self.diamond = self.add_tool("metagenomic.diamond")
                opts = {
                    "query": self.option("in_fasta"),
                    "query_type": self.option("fasta_type"),
                    "blast": self.option("blast"),
                    "database": self.option("database"),
                    "top_num": self.option("top_num"),
                    "align_len": self.option("align_len"),
                    "identity": self.option("identity"),
                    "evalue": self.option("evalue"),
                    "reference_type": "prot"
                }
                if self.option("database") in ['Custom']:
                    opts['ref'] = self.merge_fa.output_dir + "/all.scaf.fa"
                self.diamond.set_options(opts)
                self.diamond.on("end", self.set_output)
                self.diamond.run()
        else:
            if self.option("method") in ["blast"]:
                for i in os.listdir(self.split_seq.output_dir + "/" + self.option('sample') + "_fasta"):
                    n = 0
                    blast = self.add_tool("metagenomic.blast")
                    opts = {
                        "query": self.option("in_fasta"),
                        "query_type": self.query_type,
                        "blast": self.option("blast"),
                        "database": self.option("database"),
                        "reference_type": self.reference_type,
                        "top_num": self.option("top_num"),
                        "align_len": self.option("align_len"),
                        "identity": self.option("identity"),
                        "evalue": self.option("evalue"),
                    }
                    if self.option("database") in ['Custom']:
                        opts['ref'] = self.merge_fa.output_dir + "/all.scaf.fa"
                    blast.set_options(opts)
                    self.modules.append(blast)
                    n += 1
                self.logger.info(self.modules)
                if len(self.modules) > 1:
                    self.on_rely(self.modules, self.run_merge)
                elif len(self.modules) == 1:
                    self.modules[0].on("end", self.run_merge)
                for module in self.modules:
                    module.run()
                    
            elif self.option("method") in ["diamond"]:
                for i in os.listdir(self.split_seq.output_dir + "/" + self.option('sample') + "_fasta"):
                    n = 0
                    diamond = self.add_tool("metagenomic.diamond")
                    opts = {
                        "query": self.split_seq.output_dir + "/" + self.option('sample') + "_fasta/" + i ,
                        "query_type": self.option("fasta_type"),
                        "blast": self.option("blast"),
                        "database": self.option("database"),
                        "top_num": self.option("top_num"),
                        "align_len": self.option("align_len"),
                        "identity": self.option("identity"),
                        "evalue": self.option("evalue"),
                        "reference_type": "prot"
                    }
                    if self.option("database") in ['Custom']:
                        opts['ref'] =self.merge_fa.output_dir + "/all.scaf.fa"
                    diamond.set_options(opts)
                    self.modules.append(diamond)
                    n += 1
                self.logger.info(self.modules)
                if len(self.modules) > 1:
                    self.on_rely(self.modules, self.run_merge)
                elif len(self.modules) == 1:
                    self.modules[0].on("end", self.run_merge)
                for module in self.modules:
                    module.run()

    def run_merge(self):
        if os.path.exists(self.work_dir + "/result"):
            shutil.rmtree(self.work_dir + "/result")
        os.mkdir(self.work_dir + "/result")
        if len(self.modules) > 1:
            for module in self.modules:
                for i in os.listdir(module.output_dir):
                    os.link(module.output_dir + "/" + i, self.work_dir + "/result/" + i)
        elif len(self.modules) == 1:
            for i in os.listdir(self.modules[0].output_dir):
                self.logger.info(self.modules[0].output_dir + "/" + i)
                os.link(self.modules[0].output_dir + "/" + i, self.work_dir + "/result/" + i)
        add_merge(self.work_dir + "/result", "diamond", self.output_dir + "/" + self.option("sample") + ".diamond.m8.xls")
        self.end()

    def get_num(self, file):
        with open(file, "r") as f:
            lines = f.readlines()
        return len(lines)

    def run(self):
        """
        运行
        :return:
        """
        super(FaReadsModule, self).run()
        self.num = self.get_num(self.option('in_fasta').prop['path'])
        if self.num <= 500000:
            if self.option("database") in ['Custom']:
                self.merge_fa.on("end", self.run_align)
                self.run_merge_fa()
            else:
                self.run_align()         
        else:
            if self.option("database") in ['Custom']:
                self.merge_fa.on("end", self.run_align)
                self.split_seq.on("end",self.run_merge_fa)
                self.run_splitfq()
            else:
                self.split_seq.on("end", self.run_align)
                self.run_splitfq()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if self.num <= 500000:
            if self.option("method") in ["blast"]:
                link_dir(self.blast.output_dir, self.output_dir)
            elif self.option("method") in ["diamond"]:
                link_dir(self.diamond.output_dir, self.output_dir)
        else:
            pass

    def end(self):
        super(FaReadsModule, self).end()