# -*- coding: utf-8 -*-
# __author__ = 'xuting'
from __future__ import division
import math
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.files.sequence.fasta import FastaFile
from biocluster.core.exceptions import OptionError


class SamplesInfoAgent(Agent):
    """
    version 1.0
    author: xuting
    last_modify: 2015.11.06
    """
    def __init__(self, parent):
        super(SamplesInfoAgent, self).__init__(parent)
        options = [
            {"name": "fasta_path", "type": "infile", "format": "sequence.fasta_dir"}]  # 输入文件夹
        self.add_option(options)
        self.step.add_steps("sample_info_stat")
        self.on('start', self.start_info_stat)
        self.on('end', self.end_info_stat)

    def start_info_stat(self):
        self.step.sample_info_stat.start()
        self.step.update()

    def end_info_stat(self):
        self.step.sample_info_stat.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("fasta_path").is_set:
            raise OptionError("参数fasta_path不能为空")
        return True

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["samples_info", "", "样本信息结果目录"],
            ["samples_info/samples_info.txt", "xls", "样本信息统计文件"]
        ])
        super(SamplesInfoAgent, self).end()

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        total = 0
        for f in self.option("fasta_path").prop["fasta_fullname"]:
            total += os.path.getsize(f)
        total = total / (1024 * 1024 * 1024)
        total = total * 4
        total = math.ceil(total)
        self._memory = '{}G'.format(int(total))


class SamplesInfoTool(Tool):
    def __init__(self, config):
        super(SamplesInfoTool, self).__init__(config)
        self._version = 1.0

    def create_table(self):
        """
        生成samples_info表
        """
        self.logger.info('生成fasta文件夹')
        output_dir = os.path.join(self.work_dir, 'fasta')
        self.option('fasta_path').get_full_info(output_dir)
        self.logger.info('成功生成fasta文件夹,开始统计样本信息')
        sample_info_dir = os.path.join(self.work_dir, 'output/samples_info')
        if not os.path.exists(sample_info_dir):
            os.mkdir(sample_info_dir)
        file_name = os.path.join(sample_info_dir, "samples_info.txt")
        with open(file_name, "w") as f:
            head = ["sample", "reads", "bases", "avg", "min", "max"]
            f.write("\t".join(head) + "\n")
            for fasta in self.option("fasta_path").prop["fasta_fullname"]:
                fastafile = FastaFile()
                fastafile.set_path(fasta)
                if fastafile.check():
                    info_ = list()
                    fastafile.get_info()
                    basename = fastafile.prop["basename"]
                    s_name = re.sub(r"(.+)\.(fa|fasta)$", r"\1", basename)
                    info_.append(s_name)
                    info_.append(fastafile.prop["seq_number"])
                    info_.append(fastafile.prop["bases"])
                    avg = int(fastafile.prop["bases"]) / int(fastafile.prop["seq_number"])
                    avg = str(avg)
                    info_.append(avg)
                    info_.append(fastafile.prop["shortest"])
                    info_.append(fastafile.prop["longest"])
                    f.write("\t".join(info_) + "\n")
        self.logger.info('样本信息统计完毕！')

    def run(self):
        """
        运行
        """
        super(SamplesInfoTool, self).run()
        self.create_table()
        self.logger.info("退出样本信息统计模块")
        self.end()
