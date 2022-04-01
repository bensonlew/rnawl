#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == gao.hao

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class ReadsUnzipModule(Module):
    """
    解压原始数据
    """
    def __init__(self, work_id):
        super(ReadsUnzipModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},
        ]
        self.add_option(options)
        self.samplesi = ""
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
        super(ReadsUnzipModule, self).run()
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
                    elif line[2] in sample[line[1]].keys():
                        sample[line[1]][line[2]] = sample[line[1]][line[2]] + " " + sample_path
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
                if d == "r":
                    direct = "2"
                elif d == "l":
                    direct = "1"
                elif d == "s": ##add by qingchen.zhang@20200415
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
                self.linkdir(self.work_dir + "/ungiz_dir", "data")
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
                    if d ==  "s":##add by qingchen.zhang@20200415
                        direct = "s"
                    fq_name = sample + "." + direct + ".fq"
                    w.write(fq_name + '\t' + sample +'\t' + d + '\n')
        self.end()

    def linkdir(self, dirpath, dirname):
        """
		link一个文件夹下的所有文件到本module的output目录
		:param dirpath: 传入文件夹路径
		:param dirname: 新的文件夹名称
		:return:
		"""
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))
