# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.toolapps.common import link_dir


class SingleSortmernaModule(Module):
    """
    单个样品的数据处理，包含文件压缩
    """
    def __init__(self, work_id):
        super(SingleSortmernaModule, self).__init__(work_id)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "sample", "type": "string"},
            {"name": "database", "type": 'string',"default": "rfam_5.8s,rfam_5s,arc_16s,arc_23s,bac_16s,bac_23s,euk_18s,euk_28s"},
        ]
        self.add_option(options)
        self.sortmerna = self.add_tool("metagenomic.sortmerna")
        self.modules =[]

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('read1').is_set:
            raise OptionError('必须输入read1文件')
        if not self.option('read2').is_set:
            raise OptionError('必须输入read2文件')

    def run_sortmerna(self):
        opts = ({
            "read1": self.option("read1"),
            "read2": self.option("read2"),
            "sample": self.option("sample"),
            "database": self.option("database"),
        })
        self.sortmerna.set_options(opts)
        self.sortmerna.run()

    def run_gz(self):
        files = os.listdir(self.sortmerna.output_dir)
        for file in files:
            gz = self.add_tool("sequence.zip")
            opts = ({
                "file_path": self.sortmerna.output_dir + "/" + file,
                "method": "gz",
            })
            gz.set_options(opts)
            self.modules.append(gz)
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
        super(SingleSortmernaModule, self).run()
        self.sortmerna.on("end", self.run_gz)
        self.run_sortmerna()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(os.listdir(self.output_dir)) >= 1:
            for i in os.listdir(self.output_dir):
                os.remove(self.output_dir + "/" + i)
        for module in self.modules:
            link_dir(module.output_dir, self.output_dir)
        self.end()

    def end(self):
        super(SingleSortmernaModule, self).end()