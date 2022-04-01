# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re
import time, shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir


class SampleMergeBamModule(Module):
    """
    按样品合并bam文件
    author: gaohao
    last_modify: 2019.07.029
    """
    def __init__(self, work_id):
        super(SampleMergeBamModule, self).__init__(work_id)
        options = [
            {"name": "bam_dir", "type": "infile", "format": "metagbin.bam_dir"},  # bam的文件夹
            {'name': 'specimen_info', 'type': 'infile', 'format': 'meta_genomic.specimen_info'},
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
        samples = self.get_samples()
        files = os.listdir(self.option('bam_dir').prop['path'])
        n = 0
        for sample in samples:
            list = []
            for k in files:
                if k.endswith('.bam') and re.search(sample, k):
                    file = self.option("bam_dir").prop['path'] + "/" + k
                    list.append(file)
            des = " ".join(list)
            bam_sort = self.add_tool('metagbin.sample_sort_bam')
            opts = {
                "bam": des,
                "sample": sample,
                }
            bam_sort.set_options(opts)
            self.modules.append(bam_sort)
            n += 1
        self.logger.info(self.modules)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        else:
            self.modules[0].on("end", self.set_output)
        for module in self.modules:
            module.run()

    def get_samples(self):
        list = []
        with open (self.option("specimen_info").prop['path'], 'r') as f:
            lines =f.readlines()
            for line in lines[1:]:
                line =line.strip().split('\t')
                list.append(line[0])
        return list

    def run(self):
        """
        运行
        :return:
        """
        super(SampleMergeBamModule, self).run()
        self.run_bamsort()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(self.modules) > 1:
            for module in self.modules:
                link_dir(module.output_dir, self.output_dir + "/bam_sort")
        else:
            link_dir(self.modules[0].output_dir, self.output_dir + "/bam_sort")
        self.end()

    def end(self):
        super(SampleMergeBamModule, self).end()