# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171130

"""DynamicTrim.pl 用于miRNA质控去低值"""
import os
import glob
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class DynamicTrimAgent(Agent):
    def __init__(self, parent=None):
        super(DynamicTrimAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "datasplit.fastq"},  # 样本fastq文件
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 样本fastq文件夹
            {"name": "phred_score", "type": "string", "default": "20"},  # Phred得分（在0到40之间）
            {"name": "list_file", "type": "infile", "format": "sequence.file_sample"},  # fastq文件夹对应的list.txt文件
            {"name": "out_fastq", "type": "outfile", "format": "datasplit.fastq"}  # 输出结果
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("fastq").is_set and not self.option("fastq_dir").is_set:
            raise OptionError("请设置fastq文件或者fastq_dir文件夹")
        if self.option("fastq_dir").is_set and not self.option("fastq_dir").is_set:
            raise OptionError("请设置fastq_dir文件夹时设置list_file")

    def set_resource(self):
        # self._cpu = "4"
        # self._memory = "20G"
        self._cpu = "1"
        self._memory = "5G"

    def end(self):
        super(DynamicTrimAgent, self).end()


class DynamicTrimTool(Tool):
    def __init__(self, config):
        super(DynamicTrimTool, self).__init__(config)
        self._version = 1.0
        self.perl = "program/perl-5.24.0/bin/perl"
        self.dynamic_trim = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/DynamicTrim.pl"

    def run_single_dynamic_trim(self):
        """DynamicTrim.pl, 去低值"""
        trim_path = os.path.join(self.work_dir, os.path.basename(self.option("fastq").prop["path"]) + ".trimmed")
        if os.path.exists(trim_path):
            os.remove(trim_path)
        cmd = "{} {} {} -h {} -bwa".format(self.perl, self.dynamic_trim, self.option("fastq").prop["path"],\
               self.option("phred_score"))
        self.logger.info(cmd)
        command = self.add_command("dynamic_trim", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行DynamicTrim.pl完成")
        else:
            self.set_error("运行DynamicTrim.pl出错")

    def run_many_dynamic_trim(self):
        """多个fastq文件，用DynamicTrim.pl, 去低值"""
        dir_path = self.option("fastq_dir").prop["path"]
        fastqs = []
        with open(self.option("list_file").prop["path"], "r") as f:
            for line in f:
                line = line.strip().split("\t")
                if os.path.exists(line[0]):
                    fastqs.append(line[0])
                else:
                    fastqs.append(os.path.join(dir_path, line[0]))
        i = 0
        cmd_list = []
        for fastq in fastqs:
            i += 1
            cmd = "{} {} {} -h {} -bwa".format(self.perl, self.dynamic_trim, fastq, self.option("phred_score"))
            command = self.add_command("dynamic_trim_{}".format(i), cmd).run()
            cmd_list.append(command)
        self.wait()
        for command in cmd_list:
            if command.return_code == 0:
                self.logger.info("运行DynamicTrim.pl完成:{}".format(command.name))
            else:
                self.set_error("运行DynamicTrim.pl出错:{}".format(command.name))

    def set_output(self):
        files = glob.glob(r'*.trimmed')
        for f in files:
            if os.path.exists(os.path.join(self.output_dir, f)):
                os.remove(os.path.join(self.output_dir, f))
            os.link(os.path.join(self.work_dir, f), os.path.join(self.output_dir, f))
            self.option("out_fastq", os.path.join(self.output_dir, f))
        self.logger.info("Done!")

    def run(self):
        super(DynamicTrimTool, self).run()
        if self.option("fastq_dir").is_set:
            self.run_many_dynamic_trim()
        if self.option("fastq").is_set:
            self.run_single_dynamic_trim()
        self.set_output()
        self.end()
