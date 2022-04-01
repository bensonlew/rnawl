#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import subprocess
import os,re
from mbio.packages.metagenomic.common import link_file



class BacReadStatistic2Agent(Agent):
    """
    用于做微生物基因组二代数据fastq做q20q30统计
    不区分质控前和质控后
    version 1.0
    author: 顾海东
    last_modify: 20190410
    """

    def __init__(self, parent):
        super(BacReadStatistic2Agent, self).__init__(parent)
        options = [
            {"name": "q1", "type": "infile", "format": "sequence.fastq", "required": True},  # 样品左端序列
            {"name": "q2", "type": "infile", "format": "sequence.fastq", "required": True},  # 样品右端序列
            {"name": "prefix", "type": "string"}  # 输出结果名称的前缀
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '25G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./qc_stat.xls", "xls", "质控信息统计表"]
        ])
        super(BacReadStatistic2Agent, self).end()


class BacReadStatistic2Tool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(BacReadStatistic2Tool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.shell_path = "/program/sh"
        self.stat_path = self.config.PACKAGE_DIR + "/bacgenome/fastq_statisticsv2.pl"
        self.java_path = self.config.SOFTWARE_DIR +  "/program/sun_jdk1.8.0/bin/java"
        self.FastqStat_path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/QC/FastqTotalHighQualityBase.jar"
        self.bash_path = self.config.PACKAGE_DIR + "/bacgenome/q20v2.sh"
        self.bash2_path = self.config.PACKAGE_DIR + "/bacgenome/q30v2.sh"
        self.prefix = self.option("prefix") if self.option("prefix") else os.path.basename(self.option("q1").prop["path"])

    def run_q20(self):
        stat_file = self.work_dir + '/' + self.prefix + '.' + "q20.stats"
        if os.path.isfile(stat_file):
            os.remove(stat_file)
        cmd = "{} {} {} {} {} {} {}".format(self.shell_path, self.bash_path, self.java_path, self.FastqStat_path,
                                                  self.option("q1").prop["path"],self.option("q2").prop["path"], stat_file)
        self.logger.info(cmd)
        self.logger.info("开始运行q20统计")
        command = self.add_command("run_q20", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_q20完成")
        else:
            self.set_error("运行run_q20运行出错!")
            return False


    def run_q30(self):
        stat_file = self.work_dir + '/' + self.prefix + '.' + "q30.stats"
        if os.path.isfile(stat_file):
            os.remove(stat_file)
        cmd = "{} {} {} {} {} {} {}".format(self.shell_path, self.bash2_path, self.java_path, self.FastqStat_path,
                                                  self.option("q1").prop["path"],self.option("q2").prop["path"], stat_file)
        self.logger.info(cmd)
        self.logger.info("开始运行q30统计")
        path = "run_q30"
        command = self.add_command(path, cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_q30完成")
        else:
            self.set_error("运行run_q30运行出错!")
            return False

    def fastq_stat(self):
            q20 =self.work_dir + '/' + self.prefix + '.' + "q20.stats"
            q30 = self.work_dir + '/' + self.prefix + '.' + "q30.stats"
            stat_file = os.path.join(self.work_dir, "qc_stat.xls")
            cmd = "{} {} {} {} {} {} {}".format(self.perl_path, self.stat_path, self.prefix, q20, q30, self.option("q1").prop["path"], stat_file)
            self.logger.info(cmd)
            self.logger.info("开始运行汇总统计")
            path = "summary"
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行fastq_stat完成" )
            else:
                self.set_error("运行fastq_stat出错!")
                return False


    def set_output(self):
        self.logger.info("set output")
        link_file(self.work_dir+"/qc_stat.xls", self.output_dir+"/qc_stat.xls")
        self.logger.info("done")
        self.end()

    def run(self):
        """
        运行
        """
        super(BacReadStatistic2Tool, self).run()
        self.run_q20()
        self.run_q30()
        self.fastq_stat()
        self.set_output()



