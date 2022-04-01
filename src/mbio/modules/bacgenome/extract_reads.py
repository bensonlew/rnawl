# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/11'

import os
import shutil
import gevent
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_file, link_dir


class ExtractReadsModule(Module):
    """
    抽取一定乘数的reads
    """

    def __init__(self, work_id):
        super(ExtractReadsModule, self).__init__(work_id)
        option = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "fastq_list", "type": "infile", "format": "bacgenome.simple_file"},  # list表，第一列为fastq前缀，第二列为碱基数(单位bp), 第三列为基因组大小(单位兆)
            {"name": "depth_ctrl", "type": "bool", "default": True},  # 是否需要做数据抽取
            {"name": "depth_num", "type": "int", "default": 150, "min": 10}  # 抽取数据的乘数
        ]
        self.add_option(option)
        self.run_tools = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        # edit options check
        return True

    def run(self):
        super(ExtractReadsModule, self).run()
        self.run_check()

    def run_check(self):
        with open(self.option("fastq_list").prop["path"]) as file:
            lines = file.readlines()[1:]
            for line in lines:
                line = line.strip().split("\t")
                prefix, base, gsize = line
                self.logger.info(prefix)
                if self.option("depth_ctrl"):
                    self.check_then_run_extract(prefix, base, gsize)
                else:
                    self.logger.info("direct put output")
                    self.direct_put_output(prefix)
        if len(self.run_tools) > 0:
            self.on_rely(self.run_tools, self.set_output)
            for tool in self.run_tools:
                tool.run()
        else:
            gevent.spawn_later(5, self.end)

    def check_then_run_extract(self, prefix, base, gsize):
        gsize_in_bp = float(gsize) * 1204*1204
        if float(base) / gsize_in_bp > self.option("depth_num"):
            scale = self.option("depth_num") * gsize_in_bp / float(base)
            self.logger.info(scale)
            self.run_seqtk(prefix, scale)
        else:
            self.direct_put_output(prefix)

    def run_seqtk(self, prefix, scale):
        file1 = prefix + ".clean.1.fq"
        file2 = prefix + ".clean.2.fq"
        tool = self.add_tool("bacgenome.seqtk")
        tool.set_options({
            "fastq": os.path.join(self.option("fastq_dir").prop["path"], file1),
            "outfastq": file1,
            "scale": scale
        })
        self.run_tools.append(tool)
        tool2 = self.add_tool("bacgenome.seqtk")
        tool2.set_options({
            "fastq": os.path.join(self.option("fastq_dir").prop["path"], file2),
            "outfastq": file2,
            "scale": scale
        })
        self.run_tools.append(tool2)

    def direct_put_output(self, prefix):
        file1 = prefix + ".clean.1.fq"
        file2 = prefix + ".clean.2.fq"
        link_file(os.path.join(self.option("fastq_dir").prop["path"], file1), os.path.join(self.output_dir, file1))
        link_file(os.path.join(self.option("fastq_dir").prop["path"], file2), os.path.join(self.output_dir, file2))

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for tool in self.run_tools:
            link_dir(tool.output_dir, self.output_dir)
        self.logger.info("设置seqtk结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(ExtractReadsModule, self).end()
