#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class FastqStatAgent(Agent):
    """
    用于做fastq序列基本信息统计的工具
    """

    def __init__(self, parent):
        super(FastqStatAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "datasplit.fastq"},  # fastq文件
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},  # 要进行质控的文件的list,第一列文件路径，二列样本，三列r/l
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option('fastq').is_set and not self.option("list_file").is_set:
            raise OptionError("请设置fastq参数或者list_file参数")

    def set_resource(self):
        self._cpu = 10
        self._memory = '50G'

    def end(self):
        super(FastqStatAgent, self).end()


class FastqStatTool(Tool):
    def __init__(self, config):
        super(FastqStatTool, self).__init__(config)
        self.FastqStat_path = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/FastqStat.jar"
        self.java_path = self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin/java"
        self.FastqStat_sh_path = "bioinfo/seq/scripts/fastq_stat.sh"

    def fastq_stat(self):
        self.get_list_file()
        stat_file = os.path.join(self.work_dir, "fastq_stat.xls")
        cmd = "{} {} {} {} {} {}".format(self.FastqStat_sh_path, self.java_path, self.FastqStat_path, "fq_list_for_FastqStat", 10, stat_file)
        command = self.add_command("fastqstat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行FastqStat.jar完成")
            self.set_output()
        else:
            self.set_error("运行FastqStat.jar运行出错!")

    def get_list_file(self):
        if self.option("fastq").is_set:
            self.fastq_name = self.option("fastq").prop["path"].split("/")[-1]
            with open("fq_list_for_FastqStat", "wb") as w:
                w.write("{}\t{}".format(self.fastq_name, self.option("fastq").prop["path"]))
        else:
            sample_file = {}
            list_info = self.option("list_file").prop["path"]
            with open(list_info, "rb") as l, open("fq_list_for_FastqStat", "wb") as w:
                for line in l:
                    line = line.strip().split()
                    if line[1] not in sample_file:
                        sample_file[line[1]] = [line[0]]
                    else:
                        sample_file[line[1]].append(line[0])
                for s in sample_file.keys():
                    if len(sample_file[s]) == 2:
                        w.write("{}\t{}\t{}\n".format(s, sample_file[s][0], sample_file[s][1]))
                    elif len(sample_file[s]) == 1:
                        w.write("{}\t{}\n".format(s, sample_file[s][0]))

    def set_output(self):
        self.logger.info("set output")
        if os.path.exists(self.output_dir + "/fastq_stat.xls"):
            os.remove(self.output_dir + "/fastq_stat.xls")
        os.link(self.work_dir+"/fastq_stat.xls", self.output_dir+"/fastq_stat.xls")
        self.logger.info("done")
        self.end()

    def run(self):
        super(FastqStatTool, self).run()
        self.fastq_stat()
