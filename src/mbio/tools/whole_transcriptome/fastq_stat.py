# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,shicaiping'

import os
import re

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class FastqStatAgent(Agent):
    """
    用于做fastq序列基本信息统计的工具
    version 1.0
    author: qindanhua
    last modified by shicaiping at 20180507
    """

    def __init__(self, parent):
        super(FastqStatAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},
            {"name": "quality", "type": "int", "default": 33}
        ]
        self.add_option(options)
        self._memory_increase_step = 16
        self.step.add_steps('fastq_stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.fastq_stat.start()
        self.step.update()

    def step_end(self):
        self.step.fastq_stat.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("请说明序列类型，PE or SE?", code="33705201")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '32G'

    def end(self):
        super(FastqStatAgent, self).end()


class FastqStatTool(Tool):
    def __init__(self, config):
        super(FastqStatTool, self).__init__(config)
        self.FastqStat_path = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/FastqStat.jar"
        self.fastq_name = self.option("fastq").prop["path"].split("/")[-1]
        self.java_path = "program/sun_jdk1.8.0/bin/"

    def run(self):
        super(FastqStatTool, self).run()
        self.fastq_stat()

    def fastq_stat(self):
        self.get_list_file()
        cmd = "{}java -jar {} -i {} -t {} -Q {}".format(self.java_path, self.FastqStat_path, "fq_list_for_FastqStat",
                                                        10, self.option("quality"))
        self.logger.info(cmd)
        self.logger.info("开始运行FastqStat.jar")
        command = self.add_command("fastqstat", cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行FastqStat.jar完成")
            self.set_output()
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("运行FastqStat.jar运行出错!", code="33705202")

    def get_list_file(self):
        if self.option("fastq").format == "sequence.fastq":
            with open("fq_list_for_FastqStat", "wb") as w:
                w.write("{}\t{}".format(self.fastq_name, self.option("fastq").prop["path"]))
        elif self.option("fastq").format == "sequence.fastq_dir":
            sample_file = {}
            fq_dir = self.option("fastq").prop["path"]
            list_info = os.path.join(fq_dir, "list.txt")
            with open(list_info, "rb") as r, open("fq_list_for_FastqStat", "wb") as w:
                for line in r:
                    line = line.strip().split()
                    if line[1] not in sample_file:
                        sample_file[line[1]] = [line[0]]
                    else:
                        sample_file[line[1]].append(line[0])
                self.logger.info(sample_file)
                for i in sorted(sample_file):
                    if len(sample_file[i]) == 2:
                        w.write("{}\t{}\t{}\n".format(i, os.path.join(fq_dir, sample_file[i][0]),
                                                      os.path.join(fq_dir, sample_file[i][1])))
                    elif len(sample_file[i]) == 1:
                        w.write("{}\t{}\n".format(i, os.path.join(fq_dir, sample_file[i][0])))

    def set_output(self):
        self.logger.info("set output")
        os.system('rm -rf ' + self.output_dir)
        os.system('mkdir ' + self.output_dir)
        sample_list = []
        list_info = os.path.join(self.option("fastq").prop["path"], "list.txt")
        fastqstat = self.work_dir + '/fastqstat.o'
        output = self.output_dir + '/{}_fastq_stat.xls'.format(self.fastq_name)
        with open(list_info, "r") as f:
            for line in f:
                tmp = line.strip().split()
                if tmp[1] not in sample_list:
                    sample_list.append(tmp[1])
                else:
                    pass
        with open(fastqstat, "r") as r, open(output, "w") as w:
            for line in r:
                if re.search(r'^#Sample_ID', line):
                    w.write(line)
                else:
                    line1 = line.strip().split()
                    if line1[0] in sample_list:
                        w.write(line)
        self.logger.info("done")
        self.end()


