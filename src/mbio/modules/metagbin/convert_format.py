# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
import shutil
from mbio.packages.metagbin.common_function import link_dir


class ConvertFormatModule(Module):
    def __init__(self, work_id):
        super(ConvertFormatModule, self).__init__(work_id)
        options = [
            {"name": "sam", "type": "infile", "format":"metagbin.sam_dir"}, #输入比对的文件夹
            {"name": "ref_fa", "type": "infile", "format":"sequence.fasta"}, #输入的bin参考序列的文件夹
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "analysis", "type": "string", "default": "metagbin"},
        ]
        self.add_option(options)
        self.convert_format_tools = []
        self.merge = self.add_tool('metagbin.merge_bam')
        self.sortbam = self.add_module('metagbin.sort_bam')
        self.convert_result_path = ''

    def check_options(self):
        if not self.option("sam").is_set:
            raise OptionError("必须设置参数sam")
        return True

    def run_convert_format(self):
        """
        将Sam文件抽取，并对抽取的数据进行排序
        :return:
        """
        self.convert_result_path = os.path.join(self.work_dir, "Covert_format")
        if os.path.exists(self.convert_result_path):
            pass
        else:
            os.mkdir(self.convert_result_path)
        sam_dir_path = self.option('sam').prop['path']
        for sam in os.listdir(sam_dir_path):
            sam_path = os.path.join(sam_dir_path, sam)
            convert_format = self.add_tool('metagbin.convert_format')
            convert_format.set_options({
                'sam': sam_path,
                'analysis': self.option('analysis')
            })
            self.convert_format_tools.append(convert_format)
        if self.option("analysis") in ['metagbin']:
            if len(self.convert_format_tools) > 1:
                self.on_rely(self.convert_format_tools, self.run_merge)
            else:
                self.convert_format_tools[0].on('end', self.run_merge)
            for tool in self.convert_format_tools:
                tool.run()
        else:
            if len(self.convert_format_tools) > 1:
                self.on_rely(self.convert_format_tools, self.run_bam)
            else:
                self.convert_format_tools[0].on('end', self.run_bam)
            for tool in self.convert_format_tools:
                tool.run()


    def run_merge(self):
        """
        对提取结果进行merge，并转为fasta和fastq文件格式
        :return:
        """
        if self.option("analysis") in ['metagbin']:
            sort_dir = self.convert_result_path
            for i in self.convert_format_tools:
                for f in os.listdir(i.output_dir):
                    if os.path.splitext(f)[1] == '.bam':
                        file_path = os.path.join(i.output_dir, f)
                        new_path = os.path.join(self.convert_result_path, os.path.basename(file_path))
                        if os.path.exists(new_path):
                            os.remove(new_path)
                        os.link(file_path, new_path)
            self.merge.set_options({
               'sort_file': sort_dir,
               'ref_fa': self.option('ref_fa'),
               'fastq_dir': self.option('fastq_dir')
           })
            self.merge.on('end', self.set_output)
            self.merge.run()

    def run_bam(self):
        self.sort_dir = self.convert_result_path
        for i in self.convert_format_tools:
            for f in os.listdir(i.output_dir):
                if os.path.splitext(f)[1] == '.bam':
                    file_path = os.path.join(i.output_dir, f)
                    new_path = os.path.join(self.convert_result_path, os.path.basename(file_path))
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)
        self.sortbam.set_options({
            "bam_dir": self.sort_dir,
        })
        self.sortbam.on("end", self.set_output)
        self.sortbam.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
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
                oldfile_basename = os.path.basename(oldfiles[i])
                self.linkdir(oldfiles[i], os.path.join(newdir, oldfile_basename))

    def set_output(self):
        """
        设置结果文件
        :return:
        """
        self.logger.info('设置结果目录')
        if not self.option("analysis") in ['metagbin']:
            if os.path.exists(self.output_dir + '/bam_sort'):
                shutil.rmtree(self.output_dir + '/bam_sort')
            shutil.copytree(self.sortbam.output_dir + "/bam_sort", self.output_dir + '/bam_sort')
        else:
            self.linkdir(self.merge.output_dir, self.output_dir)
        self.logger.info('设置结果目录成功')
        self.end()

    def run(self):
        super(ConvertFormatModule, self).run()
        self.run_convert_format()

    def end(self):
        super(ConvertFormatModule, self).end()
