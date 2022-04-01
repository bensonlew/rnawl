#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == gao.hao

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import shutil


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
                        sample[line[1]] = {sample_path : line[2]}
                    #elif line[2] in sample[line[1]].keys():
                    #    sample[line[1]][line[2]] = sample[line[1]][line[2]] + " " + sample_path
                    else:
                        sample[line[1]][sample_path] = line[2]
                else:
                     raise OptionError('list.txt文件格式有误')
        return sample

    def run_ungiz(self):
        self.samples = self.get_list()
        samples = self.samples
        reslut_path = os.path.join(self.work_dir, "ungiz_dir")
        if os.path.exists(reslut_path):
            shutil.rmtree(reslut_path)
        os.mkdir(reslut_path)
        # 为了兼容解压之后名称相同的fq文件
        self.unzip_path = os.path.join(self.work_dir, "ungiz_dir_gz")
        if os.path.exists(self.unzip_path):
            shutil.rmtree(self.unzip_path)
        os.mkdir(self.unzip_path)
        self.fq_path = reslut_path = os.path.join(self.work_dir, "ungiz_dir_ungz")
        if os.path.exists(reslut_path):
            shutil.rmtree(reslut_path)
        os.mkdir(reslut_path)
        for sample in samples:
            for sample_path in samples[sample]:
                if sample_path.endswith(".gz"):
                    d = samples[sample][sample_path]
                    if d == "r":
                        direct = "2"
                    elif d == "l":
                        direct = "1"
                    elif d == "s":  ##add by qingchen.zhang@20200415
                        direct = "s"
                    else:
                        raise OptionError("序列的方向不对，必须为：l/r/s")
                    gunzip_fastq = self.add_tool('sequence.fastq_ungz')
                    gunzip_fastq.set_options({
                        "fastq": sample_path,
                        "sample_name": sample,
                        "direction": direct,
                        "result_path": self.unzip_path,
                        "pipeline": "metaasv",
                    })
                    self.tools.append(gunzip_fastq)
                else:
                    if os.path.exists(self.fq_path+"/"+sample):
                        os.remove(self.fq_path+"/"+sample)
                    os.link(sample_path,self.fq_path+"/"+sample)
        if len(self.tools) > 1:
            self.on_rely(self.tools, self.set_output)
        elif len(self.tools) == 0:
            self.set_output()
        else:
            self.tools[0].on('end', self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        for file in os.listdir(self.unzip_path):
            if os.path.exists(self.fq_path+"/"+file):
                os.system("cat {} {} > {}".format(self.unzip_path+"/"+file,self.fq_path+"/"+file,self.work_dir + "/ungiz_dir/" + file))
            else:
                os.link(self.unzip_path+"/"+file,self.work_dir + "/ungiz_dir/" + file)
        for file in os.listdir(self.fq_path):
            if os.path.exists(self.unzip_path+"/"+file):
                pass
            else:
                os.link(self.fq_path + "/" + file, self.work_dir + "/ungiz_dir/" + file)
        if os.path.exists(self.work_dir + "/ungiz_dir"):
            try:
                self.linkdir(self.work_dir + "/ungiz_dir", "data")
            except Exception, e:
                raise Exception('解压的结果linkdir时出错{}'.format(e))
        samples = self.samples
        list = os.path.join(self.output_dir, "data/list.txt")
        with open(list, "wb") as w:
            for sample in samples:
                #for d in samples[sample]:
                d = samples[sample][samples[sample].keys()[0]]
                direct = ""
                if d == "r":
                    direct = "2"
                if d == "l":
                    direct = "1"
                if d == "s":  ##add by qingchen.zhang@20200415
                    direct = "s"
                fq_name = sample
                w.write(fq_name + '\t' + sample + '\t' + d + '\n')
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
