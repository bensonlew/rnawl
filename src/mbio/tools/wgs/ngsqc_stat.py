# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.03

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class NgsqcStatAgent(Agent):
    """
    软件：ngsqc
    DNA产品线碱基质量分布、碱基错误率分布、fastq数据统计的工具
    """
    def __init__(self, parent):
        super(NgsqcStatAgent, self).__init__(parent)
        options = [
            {"name": "fastq_l", "type": "infile", "format": "sequence.fastq"},  # 左端fastq序列
            {"name": "fastq_r", "type": "infile", "format": "sequence.fastq"},  # 右端fastq序列
            {"name": "sample_name", "type": "string"}  #样本名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fastq_l").is_set:
            raise OptionError("请设置左端fastq序列", code="34504001")
        if not self.option("fastq_r").is_set:
            raise OptionError("请设置右端fastq序列", code="34504002")

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(NgsqcStatAgent, self).end()


class NgsqcStatTool(Tool):
    def __init__(self, config):
        super(NgsqcStatTool, self).__init__(config)
        self.ngsqc_path = "bioinfo/WGS/ngsqc"

    def run_ngsqc(self):
        """
        ngsqc: fastq数据统计软件
        """
        if self.option("sample_name"):
            self.sample_name = self.option("sample_name")
        else:
            self.sample_name = os.path.basename(self.option("fastq_l").prop["path"]).split(".fastq")[0]
        cmd = "{} -1 {} -2 {}".format(self.ngsqc_path, self.option("fastq_l").prop["path"], self.option("fastq_r").prop["path"])
        cmd += " -o {} -k {}".format(self.work_dir + "/qc", self.sample_name)
        self.logger.info("开始运行ngsqc")
        command = self.add_command("ngsqc", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ngsqc运行完成")
        else:
            self.set_error("ngsqc运行出错", code="34504001")

    def set_output(self):
        """
        对ngsqc的结果stat文件进行处理
        """
        files = os.listdir(self.work_dir + "/qc")
        if not files:
            with open(self.work_dir + "/ngsqc.o", "r") as f:
                lines = f.readlines()
                item = lines[1].strip()
                self.set_error("样本%s:%s，请检查", variables=(self.sample_name, item), code="34504004")
        for f in files:
            if os.path.exists(self.output_dir + "/" + f):
                os.remove(self.output_dir + "/" + f)
            if re.search(r".*stat", f):
                qc_stat = os.path.join(self.output_dir, self.sample_name + ".stat")
                with open(os.path.join(self.work_dir + "/qc", f), "r") as r, open(qc_stat, "a") as w:
                    w.write("#sampleID\tReads\tBases(bp)\tGC(%)\tQ30(%)\n")
                    lines = r.readlines()
                    item = lines[-3].strip().split("\t")
                    item1 = lines[-2].strip().split("\t")
                    item2 = lines[-1].strip().split("\t")
                    reads = int(item1[2]) + int(item2[2])
                    w.write(item[1] + "\t" + str(reads) + "\t" + item[3] + "\t" + item[9] + "\t" + item[10] + "\n")
            else:
                os.link(self.work_dir + "/qc/" + f, self.output_dir + "/" + f)

    def run(self):
        super(NgsqcStatTool, self).run()
        self.run_ngsqc()
        self.set_output()
        self.end()
