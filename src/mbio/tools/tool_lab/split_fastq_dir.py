# -*- coding: utf-8 -*-
#__author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import shutil
import os

class SplitFastqDirAgent(Agent):
    def __init__(self, parent):
        super(SplitFastqDirAgent, self).__init__(parent)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "lines", "type": "int", "default": 20000000},  # 序列数
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("in_fastq").is_set:
            raise OptionError("请传入in_fastq文件！", code="")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(SplitFastqDirAgent, self).end()


class SplitFastqDirTool(Tool):
    """
    version 1.0对fastq序列进行拆分
    """
    def __init__(self, config):
        super(SplitFastqDirTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + "/metagbin/"
        self.sample_path = {}

    def split_fastq(self):
        """
        切分fastq文件
        :return:
        """
        lines = self.option('lines')
        n = 1
        for sample in self.sample_path:
            fastq1 = self.sample_path[sample][0]
            fastq2 = self.sample_path[sample][1]
            if os.path.exists(self.work_dir + "/" + sample):
                shutil.rmtree(self.work_dir + "/" + sample)
            os.mkdir(self.work_dir + "/" + sample)
            cmd1 = "{}split_reads.sh {} {} {}".format(self.sh_path, lines, fastq1, self.work_dir + "/" + sample + "/" + sample + "_l_fastq_")
            cmd2 = "{}split_reads.sh {} {} {}".format(self.sh_path, lines, fastq2, self.work_dir + "/" + sample + "/" + sample + "_r_fastq_")
            command1 = self.add_command("run_split_fastq_{}_l".format(str(n)), cmd1).run()
            n+=1
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("设置切分结果目录成功!")
            else:
                self.set_error("设置切分结果目录运行出错!")
            command2 = self.add_command("run_split_fastq_{}_r".format(str(n)), cmd2).run()
            n+=1
            self.wait(command2)
            if command2.return_code == 0:
                self.logger.info("设置切分结果目录成功!")
            else:
                self.set_error("设置切分结果目录运行出错!")

    def set_output(self):
        with open(self.work_dir+"/sample_split_file.txt","w") as t1:
            for sample in self.sample_path:
                sample_all = {}
                for file in os.listdir(self.work_dir + "/" + sample):
                    number = file.split("_fastq_")[-1]
                    if sample + "_" + number in sample_all:
                        sample_all[sample + "_" + number].append(file)
                    else:
                        sample_all[sample + "_" + number] = [file]
                        t1.write(sample+"\t"+sample + "_" + number+"\n")
                for sample_split in sample_all:
                    if len(sample_all[sample_split]) != 2:
                        raise OptionError("fastq文件拆分后数量不一致！", code="")
                    if os.path.exists(self.output_dir + "/" + sample_split):
                        shutil.rmtree(self.output_dir + "/" + sample_split)
                    os.mkdir(self.output_dir + "/" + sample_split)
                    os.link(self.work_dir + "/" + sample + "/" + sample_all[sample_split][0],
                            self.output_dir + "/" + sample_split + "/" + sample_all[sample_split][0])
                    os.link(self.work_dir + "/" + sample + "/" + sample_all[sample_split][1],
                            self.output_dir + "/" + sample_split + "/" + sample_all[sample_split][1])
                    with open(self.output_dir + "/" + sample_split + "/list.txt", "w") as t:
                        if "_l_fastq_" in sample_all[sample_split][0]:
                            style1 = "l"
                        else:
                            style1 = "r"
                        if "_r_fastq_" in sample_all[sample_split][1]:
                            style2 = "r"
                        else:
                            style2 = "l"
                        t.write(sample_all[sample_split][0] + "\t" + sample_split + "\t" + style1 + "\n" +
                                sample_all[sample_split][1] + "\t" + sample_split + "\t" + style2 + "\n")
        self.end()

    def get_info(self):
        with open(self.option('in_fastq').prop['path'] + "/list.txt") as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                if tmp[1] in self.sample_path.keys():
                    if tmp[2] == 'l':
                        self.sample_path[tmp[1]].insert(0, self.option('in_fastq').prop['path']+ '/' + tmp[0])
                    else:
                        self.sample_path[tmp[1]].append(self.option('in_fastq').prop['path'] + '/' + tmp[0])
                else:
                    self.sample_path[tmp[1]] = [self.option('in_fastq').prop['path'] + '/' + tmp[0]]

    def run(self):
        super(SplitFastqDirTool, self).run()
        self.get_info()
        self.split_fastq()
        self.set_output()
