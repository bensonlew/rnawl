#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class FastqStatAgent(Agent):
    """
    用于做fastq序列基本信息统计的工具
    version 1.0
    author: qindanhua
    last_modify: 2016.06.23
    """

    def __init__(self, parent):
        super(FastqStatAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq,sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"}
        ]
        self.add_option(options)
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
            raise OptionError("请说明序列类型，PE or SE?")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./fastq_stat.xls", "xls", "fastq信息统计表"]
        ])
        super(FastqStatAgent, self).end()


class FastqStatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(FastqStatTool, self).__init__(config)
        self.FastqStat_path = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/FastqStat.jar"
        self.fastq_name = self.option("fastq").prop["path"].split("/")[-1]
        self.java_path = "program/sun_jdk1.8.0/bin/"

    def fastq_stat(self):
        self.get_list_file()
        cmd = "{}java -jar {} -i {} -t {}".format(self.java_path, self.FastqStat_path, "fq_list_for_FastqStat", 10)
        self.logger.info(cmd)
        self.logger.info("开始运行FastqStat.jar")
        command = self.add_command("fastqstat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行FastqStat.jar完成")
            self.set_output()
        else:
            self.logger.info("第一次运行失败，开始重新运行")
            command.rerun()
            self.wait(command)
            if command.return_code == 0:
                self.set_error("运行FastqStat.jar运行出错!")
                return False

    def get_list_file(self):
        if self.option("fastq").format == "sequence.fastq":
            with open("fq_list_for_FastqStat", "wb") as w:
                w.write("{}\t{}".format(self.fastq_name, self.option("fastq").prop["path"]))
        elif self.option("fastq").format == "sequence.fastq_dir":
            sample_file = {}
            fq_dir = self.option("fastq").prop["path"]
            list_info = os.path.join(fq_dir, "list.txt")
            with open(list_info, "rb") as l, open("fq_list_for_FastqStat", "wb") as w:
                for line in l:
                    line = line.strip().split()
                    if line[1] not in sample_file:
                        sample_file[line[1]] = [line[0]]
                    else:
                        sample_file[line[1]].append(line[0])
                self.logger.info(sample_file)
                for i in sample_file:
                    # print os.path.join(fq_dir, sample_file[i][0])
                    if len(sample_file[i]) == 2:
                        w.write("{}\t{}\t{}\n".format(i, os.path.join(fq_dir, sample_file[i][0]), os.path.join(fq_dir, sample_file[i][1])))
                    elif len(sample_file[i]) == 1:
                        w.write("{}\t{}\n".format(i, os.path.join(fq_dir, sample_file[i][0])))

    def set_output(self):
        self.logger.info("set output")
        os.system('rm -rf '+self.output_dir)
        os.system('mkdir '+self.output_dir)
        os.link(self.work_dir+'/fastqstat.o', self.output_dir+'/{}_fastq_stat.xls'.format(self.fastq_name))
        os.system("sed -i '1d' {}".format(self.output_dir+'/{}_fastq_stat.xls'.format(self.fastq_name)))
        os.system("sed -i '$d' {}".format(self.output_dir+'/{}_fastq_stat.xls'.format(self.fastq_name)))
        os.system("sed -i '$d' {}".format(self.output_dir+'/{}_fastq_stat.xls'.format(self.fastq_name)))
        self.logger.info("done")
        self.end()

    def run(self):
        """
        运行
        """
        super(FastqStatTool, self).run()
        self.fastq_stat()
        # self.set_output()
