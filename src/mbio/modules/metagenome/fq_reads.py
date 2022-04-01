# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import link_dir,add_merge


class FqReadsModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(FqReadsModule, self).__init__(work_id)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq"},
            {"name": "read_type", "type": "int", "default": 1},
            {"name": "sample", "type": "string"},
            {"name": "method", "type": "string", "default": "blast"},  # blast或者diamond
            {'name': 'blast', 'type': 'string', 'default': 'blastn'},  # blastn,blastp,blastx
            {'name': 'database', 'type': 'string', 'default': 'NT'},  # NT,NR,Swiss-Prot,Custom
            {'name': 'ref_dir', 'type': 'infile', 'format': 'sequence.fasta_dir'},  # Custom时输入文件
            {'name': 'top_num', 'type': 'string', 'default': '1'},
            {'name': 'align_len', 'type': 'string', 'default': '50'},
            {'name': 'identity', 'type': 'string', 'default': '30'},
            {'name': 'evalue', 'type': 'string', 'default': '1e-5'}
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
        if not self.option('in_fastq').is_set:
            raise OptionError('必须输入in_fastq文件')

    def run_splitfq(self):
        opts = {
            "in_fastq": self.option('in_fastq'),
            "sample": self.option('sample'),
            "read_type": self.option('read_type'),
        }
        self.split_seq.set_options(opts)
        self.split_seq.on("end",self.run_align)
        self.split_seq.run()

    def run_align(self):
        for i in os.listdir(self.split_seq.output_dir + "/" + self.option('sample') + "_fastq_" + str(self.option("read_type"))):
            n = 0
            split_seq = self.add_module("metagenome.reads_align")
            opts = {
                "method": self.option('method'),
                "reads_fq": self.split_seq.output_dir + "/" + self.option('sample') + "_fastq_" + str(self.option("read_type")) + "/" + i,
                "blast": self.option("blast"),
                "database": self.option("database"),
                "ref_dir": self.option("ref_dir"),
                "top_num": self.option("top_num"),
                "align_len": self.option("align_len"),
                "identity": self.option("identity"),
                "evalue": self.option("evalue"),
            }
            split_seq.set_options(opts)
            self.modules.append(split_seq)
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
        add_merge(self.work_dir + "/result", "", self.output_dir + "/" + self.option("sample") + "_" + str(self.option("read_type")) + "_m8.xls")
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        super(FqReadsModule, self).run()
        self.run_splitfq()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        pass

    def end(self):
        super(FqReadsModule, self).end()