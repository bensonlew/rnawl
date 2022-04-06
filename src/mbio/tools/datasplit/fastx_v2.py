# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171121

"""fastxtoolkit  用于统计碱基质量信息"""
import os
import re
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class FastxV2Agent(Agent):
    """
    统计碱基质量信息及Q20、Q30
    """
    def __init__(self, parent=None):
        super(FastxV2Agent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "datasplit.fastq"},  # fastq文件
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},  # 要进行质控的文件的list,第一列文件路径，二列样本，三列r/l
        ]
        self.add_option(options)

    def check_option(self):
        """
        参数检测
        """
        if not self.option('fastq').is_set and not self.option("list_file").is_set:
            raise OptionError("请设置fastq参数或者list_file参数")

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 5
        self._memory = '30G'


class FastxV2Tool(Tool):
    """
    """
    def __init__(self, config):
        super(FastxV2Tool, self).__init__(config)
        self._version = 1.0
        self.fastx_dir = "bioinfo/seq/fastx_toolkit_0.0.14/"
        self.python_dir = "miniconda2/bin/python"
        self.q20q30_stat = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/q20q30_stat.py"
        self.fastqs, fastqs = list(), list()
        if self.option("list_file").is_set:
            fastqs = self.option("list_file").prop["fastqs"]
        else:
            fastqs.append(self.option("fastq").prop["path"])
        for f in fastqs:
            m = re.match(r"(.+).gz|.gzip", os.path.basename(f))
            if m:
                fastq = os.path.join(self.work_dir, m.group(1))
                os.system("gunzip -c {} > {}".format(f, fastq))
                self.fastqs.append(fastq)
            else:
                self.fastqs.append(f)
        self.fastx = list()

    def fastxtoolkit(self):
        """
        统计碱基质量
        """
        cmd_list = list()
        i = 0
        log_path = os.path.join(self.work_dir, "log.sh")
        with open(log_path, 'wb') as w:
            my_str = "fastxtoolkit"
            w.write(my_str.center(79, "#"))
            w.write("\n")
        for fastq in self.fastqs:
            i += 1
            file_name = os.path.join(self.output_dir, os.path.basename(fastq) + ".fastxstat")
            self.fastx.append(file_name)
            cmd = self.fastx_dir + "fastx_quality_stats" + " -i " + fastq + " -o " + file_name
            self.logger.debug(cmd)
            with open(log_path, 'wb') as w:
                w.write(cmd)
                w.write(cmd)
            # command = self.add_command("fastx_quality_stats" + str(i), cmd)
            command = self.add_command("fastx_quality_stats" + str(i), cmd, False, 0, True)  # 非正常结束时忽略不报错
            cmd_list.append(command)
        for mycmd in cmd_list:
            self.logger.info("开始运行fastx_quality_stats")
            mycmd.run()
        self.wait()
        with open(self.work_dir + "/fastx_status.txt", "w") as w:
            self.logger.info(self.work_dir + "/fastx_status.txt")
            for mycmd in cmd_list:
                if mycmd.return_code == 0:
                    self.logger.info("运行fastx_quality_stats完成")
                else:
                    w.write("{}:运行fastx_quality_stats出错，请检查".format(mycmd.name))
                    self.logger.info("运行fastx_quality_stats出错")
                    # self.set_error("运行fastx_quality_stats出错")

    def q20_q30(self):
        cmd_list = list()
        i = 0
        for fastq in self.fastqs:
            i += 1
            file_name = os.path.join(self.output_dir, os.path.basename(fastq) + ".q20q30")
            cmd = "{} {} -i {} -o {}".format(self.python_dir, self.q20q30_stat, fastq, file_name)
            command = self.add_command("q20q30_stat" + str(i), cmd)
            cmd_list.append(command)
        for mycmd in cmd_list:
            self.logger.info("开始运行" + mycmd.name)
            mycmd.run()
        self.wait()
        for mycmd in cmd_list:
            if mycmd.return_code == 0:
                self.logger.info(mycmd.name + " 统计完成")
            else:
                self.set_error(mycmd.name + " 统计出错")

    def run(self):
        super(FastxV2Tool, self).run()
        self.fastxtoolkit()
        self.q20_q30()
        self.end()
