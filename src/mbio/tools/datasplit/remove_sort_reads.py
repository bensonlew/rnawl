#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import glob


class RemoveSortReadsAgent(Agent):
    """
    宏基因删除短的reads
    """
    def __init__(self, parent):
        super(RemoveSortReadsAgent, self).__init__(parent)
        options = [
            {"name": "fastq_r", "type": "infile", "format": "sequence.fastq"},  # 输入文件PE的右端序列
            {"name": "fastq_l", "type": "infile", "format": "sequence.fastq"},  # PE的左端序列
            {"name": "min_length", "type": "string", "default": "50"},  # 删除短于此值的reads
            {"name": "sample_name", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("fastq_r").is_set:
            raise OptionError("请传入fastq_r序列文件")
        if not self.option("fastq_l").is_set:
            raise OptionError("请传入fastq_l序列文件")
        if not self.option("sample_name"):
            raise OptionError("请设置参数sample_name")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(RemoveSortReadsAgent, self).end()


class RemoveSortReadsTool(Tool):
    def __init__(self, config):
        super(RemoveSortReadsTool, self).__init__(config)
        self.perl = "program/perl-5.24.0/bin/perl"
        self.remove_reads_path = self.config.PACKAGE_DIR + "/datasplit/remove_short_reads.pair.2.pl"
        self.pigz_path = "bioinfo/seq/pigz-2.4/pigz"

    def remove_short_reads(self):
        """
        宏基因组用脚本remove_short_reads.pair.2.pl 删除部分reads
        """
        cmd = "{} {} {}".format(self.perl, self.remove_reads_path, self.option("fastq_l").prop["path"])
        cmd += " {} {} {}".format(self.option("fastq_r").prop["path"], self.option("min_length"), self.output_dir + "/" + self.option("sample_name"))
        command = self.add_command("remove_short_reads", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("remove_short_reads运行完成")
        else:
            self.set_error("remove_short_reads运行失败")

    def file_gz(self):
        """
        对文件进行压缩
        """
        file_list = []
        for f in os.listdir(self.output_dir):
            file_list.append(os.path.join(self.output_dir, f))
        cmd = "{} -k {}".format(self.pigz_path, " ".join(file_list))
        command = self.add_command("pigz_file", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("pigz_file压缩文件运行完成")
        else:
            self.set_error("pigz_file压缩文件运行失败")

    def run(self):
        super(RemoveSortReadsTool, self).run()
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        self.remove_short_reads()
        self.file_gz()
        self.end()
