# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import get_info,link_dir


class CentrifugeModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(CentrifugeModule, self).__init__(work_id)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},#指控后的reads序列的目录
            {"name": "readlen", "type": "int", "default": 100},
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

    def run_kraken(self):
        self.sample_path = get_info(self.option('fa_dir').prop['path'])
        n=0
        for sample in sorted(self.sample_path.keys()):
            kraken = self.add_tool("metagenomic.centrifuge")
            opts = {
                "read1": self.sample_path[sample][0],
                "read2": self.sample_path[sample][1],
                "sample": sample,
                "readlen": self.option("readlen"),
                }
            kraken.set_options(opts)
            self.modules.append(kraken)
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
        super(CentrifugeModule, self).run()
        self.run_kraken()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(os.listdir(self.output_dir)) >= 1:
            for i in os.listdir(self.output_dir):
                os.remove(self.output_dir + "/" + i )
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
        super(CentrifugeModule, self).end()