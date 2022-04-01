# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
## @20200226

from biocluster.core.exceptions import OptionError
from mbio.packages.toolapps.common import link_dir
from biocluster.module import Module
import os
import shutil
import gevent



class TrinityModule(Module):
    """
    小工具：宏转录组Trinity软件组装拼接
    """
    def __init__(self, work_id):
        super(TrinityModule, self).__init__(work_id)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq序列文件夹
            {'name': 'min_contig_length', 'type': 'int', 'default': 200}, #trnity报告出的最短的Contig长度，默认为200
            {'name': 'kmer_size', "type": 'int', 'default': 25}, # Trinity软件使用的kmer的长度，默认为25
            {'name': 'min_kmer_cov', 'type': 'int', 'default': 2},  # 最小kmer的计数
        ]
        self.add_option(options)
        self.step.add_steps('sequence', 'trinity')
        self.sequence = self.add_module('metagenome.reads_unzip')
        self.tool_list = []

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('in_fastq').is_set:
            raise OptionError('请输入原始序列文件夹或者质控优质序列文件夹！', code="24400201")
        return True

    def run_sequence(self):
        """
        对数据进行解压，上传序列为压缩格式
        :return:
        """
        self.logger.info("对上传fastq文件夹进行压缩")
        opts = {
            'fastq_dir': self.option('in_fastq'),
        }
        self.sequence.set_options(opts)
        self.sequence.run()

    def run_trinity(self):
        """
        并行运行Trinity软件进行组装拼接
        :return:
        """
        self.logger.info("start trinity predict")
        samples = self.get_list()
        self.logger.info("samples:{}".format(samples))
        for sample in samples.keys():
            self.trinity = self.add_tool('assemble.trinity_assemble')
            sample_direct = samples[sample]
            if "l" in sample_direct.keys() or "r" in sample_direct.keys():
                opts = {
                    'fq_type': 'PE',
                    'fq_l': sample_direct['l'],
                    'fq_r': sample_direct['r'],
                    'sample': sample,
                    'min_contig_length': self.option("min_contig_length"),
                    'kmer_size': self.option('kmer_size'),
                    'min_kmer_cov': self.option('min_kmer_cov')
                }
                self.trinity.set_options(opts)
            elif 's' in sample_direct.keys():
                opts = {
                    'fq_type': 'SE',
                    'fq_s': sample_direct['s'],
                    'sample': sample,
                    'min_contig_length': self.option("min_contig_length"),
                    'kmer_size': self.option('kmer_size'),
                    'min_kmer_cov': self.option('min_kmer_cov')
                }
                self.trinity.set_options(opts)
            self.tool_list.append(self.trinity)
        if len(self.tool_list) >1 :
            self.on_rely(self.tool_list, self.set_output)
        else:
            self.tool_list[0].on("end", self.set_output)
        for tool in self.tool_list:
            tool.run()
            gevent.sleep(0)

    def get_list(self):
        """
        根据list.txt文件获取样品和样本名称的对应关系
        :return:
        """
        list_path = os.path.join(self.sequence.output_dir,"data","list.txt")
        if os.path.exists(list_path):
            self.logger.info(list_path)
        else:
            self.set_error("list_path不正确！", code="24400201")
        sample = {}
        with open(list_path, "rb") as l:
            for line in l:
                line = line.strip().split("\t")
                if len(line) == 3:
                    sample_path = os.path.join(self.sequence.output_dir,"data", line[0])
                    if line[1] not in sample:
                        sample[line[1]] = {line[2]: sample_path}
                    elif line[2] in sample[line[1]].keys():
                        sample[line[1]][line[2]] = sample[line[1]][line[2]] + " " + sample_path
                    else:
                        sample[line[1]][line[2]] = sample_path
                else:
                    raise OptionError('list.txt文件格式有误', code="24400202")
        return sample

    def run(self):
        """
        运行
        :return:
        """
        super(TrinityModule, self).run()
        self.sequence.on("end", self.run_trinity)
        self.run_sequence()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("开始将结果文件复制到结果目录")
        trinity = os.path.join(self.output_dir, "trinity")
        if os.path.exists(trinity):
            shutil.rmtree(trinity)
        else:
            os.mkdir(trinity)
        for tool in self.tool_list:
            link_dir(tool.output_dir, trinity)
        self.logger.info("链接结果文件完成")
        self.end()

    def end(self):
        """
        运行结束
        :return:
        """
        result_dir = self.add_upload_dir(os.path.join(self.output_dir, 'trinity'))
        self.logger.info("开始上传结果文件")
        result_dir.add_relpath_rules([
            [".", "", "Trinity组装拼接结果文件目录"],
        ])
        super(TrinityModule, self).end()