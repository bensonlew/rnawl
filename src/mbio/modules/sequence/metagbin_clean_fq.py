#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == gaohao

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class MetagbinCleanFqModule(Module):
    """
    对clean data数据解压，并生成
    """
    def __init__(self, work_id):
        super(MetagbinCleanFqModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},
        ]
        self.add_option(options)
        self.tools = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("必须输入样本文件夹！")
        else:
            return True

    def run(self):
        super(MetagbinCleanFqModule, self).run()
        self.run_ungiz()

    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop['path'], "list.txt")
        if os.path.exists(list_path):
            self.logger.info(list_path)
        sample = {}
        with open(list_path, "rb") as l:
            for line in l:
                line = line.strip().split()
                if len(line) == 3:
                    sample_path = os.path.join(self.option("fastq_dir").prop['path'], line[0])
                    if line[1] not in sample:
                        sample[line[1]] = {line[2]: sample_path}
                    else:
                        sample[line[1]][line[2]] = sample_path
                else:
                     raise OptionError('list.txt文件格式有误')
        return sample

    def run_ungiz(self):
        self.samples = self.get_list()
        samples = self.samples
        reslut_path = os.path.join(self.work_dir, "ungiz_dir")
        if not os.path.exists(reslut_path):
            os.mkdir(reslut_path)
        for sample in samples:
            for d in samples[sample]:
                direct = ""
                if d == "r":
                    direct = "2"
                elif d == "l":
                    direct = "1"
                elif d == 's':
                    direct = "s"
                else:
                    raise OptionError("序列的方向不对，必须为：l/r/s")
                gunzip_fastq = self.add_tool('sequence.fastq_ungz')
                gunzip_fastq.set_options({
                    "fastq": samples[sample][d],
                    "sample_name": sample,
                    "direction": direct,
                    "result_path": reslut_path
                })
                self.tools.append(gunzip_fastq)
        if len(self.tools) > 1:
            self.on_rely(self.tools, self.set_output)
        else:
            self.tools[0].on('end', self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        if os.path.exists(self.work_dir + "/ungiz_dir"):
            try:
                link_dir(self.work_dir + "/ungiz_dir", self.output_dir + "/data")
            except Exception, e:
                raise Exception('解压的结果linkdir时出错{}'.format(e))
        samples = self.samples
        list = os.path.join(self.output_dir, "data/list.txt")
        with open(list, "wb") as w:
            for sample in samples:
                for d in samples[sample]:
                    direct = ""
                    if d == "r":
                        direct = "2"
                    if d ==  "l":
                        direct = "1"
                    elif d == 's':
                        direct = "s"
                    fq_name = sample + "." + direct + ".fq"
                    w.write(fq_name + '\t' + sample +'\t' + d + '\n')
        self.end()
