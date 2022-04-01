# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import get_info,link_dir


class ReadsAlignDatabaseModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(ReadsAlignDatabaseModule, self).__init__(work_id)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fastq_dir"},#指控后的reads序列的目录
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
        if not self.option('fa_dir').is_set:
            raise OptionError('必须输入fa_dir文件夹')

    def run_align(self):
        self.sample_path = get_info(self.option('fa_dir').prop['path'])
        n=0
        for sample in sorted(self.sample_path.keys()):
            sample_fq = self.add_module("metagenome.sample_fq")
            opts = {
                "read1": self.sample_path[sample][0],
                "read2": self.sample_path[sample][1],
                "sample": sample,
                "method": self.option("method"),
                "blast": self.option("blast"),
                "database": self.option("database"),
                "ref_dir": self.option("ref_dir"),
                "top_num": self.option("top_num"),
                "align_len": self.option("align_len"),
                "identity": self.option("identity"),
                "evalue": self.option("evalue"),
                }
            sample_fq.set_options(opts)
            self.modules.append(sample_fq)
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
        super(ReadsAlignDatabaseModule, self).run()
        self.run_align()

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
        super(ReadsAlignDatabaseModule, self).end()