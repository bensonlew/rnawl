# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import get_info,link_dir


class SampleFqModule(Module):
    """
    按样品进行分析，按样品结果整理read1和read2的结果
    """
    def __init__(self, work_id):
        super(SampleFqModule, self).__init__(work_id)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},#
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},  #
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
        self.sample_path = defaultdict(list)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('read1').is_set:
            raise OptionError('必须输入read1文件夹')
        if not self.option('read2').is_set:
            raise OptionError('必须输入read2文件夹')

    def run_seq(self):
        n =1
        for i in [self.option("read1").prop['path'], self.option("read2").prop['path']]:
            fq_read = self.add_module("metagenome.fq_reads")
            opts = {
                "in_fastq": i,
                "sample": self.option("sample"),
                "read_type": n,
                "method": self.option("method"),
                "blast": self.option("blast"),
                "database": self.option("database"),
                "ref_dir": self.option("ref_dir"),
                "top_num": self.option("top_num"),
                "align_len": self.option("align_len"),
                "identity": self.option("identity"),
                "evalue": self.option("evalue"),
                }
            fq_read.set_options(opts)
            self.modules.append(fq_read)
            n += 1
        self.logger.info(self.modules)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(SampleFqModule, self).run()
        self.run_seq()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(os.listdir(self.output_dir)) >= 1:
            for i in os.listdir(self.output_dir):
                os.remove(self.output_dir + "/" + i)
        if len(self.modules) > 1:
            for module in self.modules:
                for i in os.listdir(module.output_dir):
                    os.link(module.output_dir + "/" + i, self.output_dir + "/" + i)
        elif len(self.modules) == 1:
            for i in os.listdir(self.modules[0].output_dir):
                self.logger.info(self.modules[0].output_dir + "/" + i)
                os.link(self.modules[0].output_dir + "/" + i, self.output_dir + "/" + i)
        self.end()

    def end(self):
        super(SampleFqModule, self).end()