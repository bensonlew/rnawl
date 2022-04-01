# -*- coding: utf-8 -*-

import os
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import get_info


class Metaphlan3Module(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(Metaphlan3Module, self).__init__(work_id)
        options = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "min_cu_len", "type": "int", "default": 2000},
            {"name": "stat", "type": "string", "default": "tavg_g"},
            {"name": "stat_q", "type": "float", "default": 0.2},
        ]
        self.modules = []
        self.add_option(options)
        self.sample_path = defaultdict(list)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fq_dir').is_set:
            raise OptionError('必须输入fa_dir文件夹')

    def run_metaphlan(self):
        self.sample_path = get_info(self.option('fq_dir').prop['path'])
        for sample in sorted(self.sample_path.keys()):
            metaphlan = self.add_tool("metagenomic.metaphlan3")
            opts = {
                "read1": self.sample_path[sample][0],
                "read2": self.sample_path[sample][1],
                "sample": sample,
                "min_cu_len": self.option("min_cu_len"),
                "stat": self.option("stat"),
                "stat_q": self.option("stat_q"),
                }
            metaphlan.set_options(opts)
            self.modules.append(metaphlan)
        self.logger.info(self.modules)
        self.on_rely(self.modules, self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(Metaphlan3Module, self).run()
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
            for i in os.listdir(self.modules[0].output_dir):
                self.logger.info(self.modules[0].output_dir + "/" + i)
                os.link(self.modules[0].output_dir + "/" + i, self.output_dir + "/" + i)
        self.end()

    def end(self):
        super(Metaphlan3Module, self).end()
