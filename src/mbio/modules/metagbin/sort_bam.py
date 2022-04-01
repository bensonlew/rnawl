# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re
import time, shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir


class SortBamModule(Module):
    """
    对多个bam文件sort、index
    author: gaohao
    last_modify: 2019.07.029
    """

    def __init__(self, work_id):
        super(SortBamModule, self).__init__(work_id)
        options = [
            {"name": "bam_dir", "type": "infile", "format": "metagbin.bam_dir"},  # bam的文件夹
        ]
        self.modules = []
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('bam_dir').is_set:
            raise OptionError('必须输入bam_dir文件夹')

    def run_bamsort(self):
        files = os.listdir(self.option('bam_dir').prop['path'])
        n = 0
        for k in files:
            bam_sort = self.add_tool('metagbin.bam_sort')
            file = self.option("bam_dir").prop['path'] + "/" + k
            opts = {
                "bam": file,
            }
            bam_sort.set_options(opts)
            self.modules.append(bam_sort)
            n += 1
        self.logger.info(self.modules)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(SortBamModule, self).run()
        self.run_bamsort()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for module in self.modules:
            link_dir(module.output_dir, self.output_dir + "/bam_sort")
        self.end()

    def end(self):
        super(SortBamModule, self).end()