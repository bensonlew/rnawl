# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import get_info,link_dir


class MetaphlanModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(MetaphlanModule, self).__init__(work_id)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},#指控后的reads序列的目录
            {"name": "read_min_len", "type": "int", "default": 70},
            {"name": "bt2_ps", "type": "string", "default": "very-sensitive"},
            # sensitive,very-sensitive,sensitive-local,very-sensitive-local
            {"name": "tax_lev", "type": "string", "default": "a"},  # a,k,p,c,o,f,g,s
            {"name": "min_cu_len", "type": "int", "default": 2000},
            {"name": "stat", "type": "string", "default": "avg_g"},  # avg_g,avg_l,tavg_g,tavg_l,wavg_g,wavg_l,med
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

    def run_metaphlan(self):
        self.sample_path = get_info(self.option('fa_dir').prop['path'])
        n=0
        for sample in sorted(self.sample_path.keys()):
            self.metaphlan = self.add_tool("metagenomic.metaphlan")
            opts = {
                "read1": self.sample_path[sample][0],
                "read2": self.sample_path[sample][1],
                "sample": sample,
                "read_min_len": self.option("read_min_len"),
                "bt2_ps": self.option("bt2_ps"),
                "tax_lev": self.option("tax_lev"),
                "min_cu_len": self.option("min_cu_len"),
                "stat": self.option("stat"),
                }
            self.metaphlan.set_options(opts)
            self.modules.append(self.metaphlan)
            self.metaphlan.on('end', 'metaphlan{}'.format(n))
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
        super(MetaphlanModule, self).run()
        self.run_metaphlan()

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
            for i in  os.listdir(self.modules[0].output_dir):
                os.link(self.modules[0].output_dir + "/" + i, self.output_dir + "/" + i)
        self.end()

    def end(self):
        super(MetaphlanModule, self).end()