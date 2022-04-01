# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20181217

import os
import re
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class UniformAgent(Agent):
    """
    fastq均一化
    """
    def __init__(self, parent=None):
        super(UniformAgent, self).__init__(parent)
        options = [
            {"name": "sample_list", "type": "infile", "format": "noref_wgs.list_file", "required": True},
            # fastq路径list.txt文件，第一列分析样本名，第二列批次样本名，第三列fastq_l路径，第四列fastq_r路径
            {"name": "analysis_method", "type": "string", "default": "ipyrad"},  # 分析方法,ipyrad/stacks
            {"name": "enzyme_method", "type": "string", "required": True},  # 酶切方案,GBS/RAD
            {"name": "uniform_length", "type": "int", "default": 120},  # -l, ngsuniform的均一长度
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("analysis_method") not in ["ipyrad", "stacks"]:
            raise OptionError("analysis_method: %s只能是ipyrad/stacks" , variables=( self.option("analysis_method")), code="35501205")
        if self.option("enzyme_method") not in ["RAD", "GBS"]:
            raise OptionError("酶切方案: %s只能是RAD/GBS" , variables=( self.option("enzyme_method")), code="35501206")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(UniformAgent, self).end()


class UniformTool(Tool):
    def __init__(self, config):
        super(UniformTool, self).__init__(config)
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.ngsuniform = self.config.SOFTWARE_DIR + "/bioinfo/noRefWGS/ngsuniform"

    def get_sample_info(self):
        """
        得到同一个样本的R1端和R2端fastq
        """
        self.sample_info = {}
        self.sample_list = []
        with open(self.option("sample_list").prop["path"], "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if item[0] not in self.sample_list:
                    self.sample_list.append(item[0])
                    self.sample_info[item[0]] = {"r": [], "l": []}
                self.sample_info[item[0]]["l"].append(item[2])
                self.sample_info[item[0]]["r"].append(item[3])

    def get_ipyrad_fastq(self):
        """
        RAD: 将同一个样本R1端的fastq cat到一起,命名方式_R1_
        GBS: 将同一个样本R1、R2端的fastq 分别cat到一起，命名方式_R1_、_R2_
        """
        for sample in self.sample_list:
            fastq_l = os.path.join(self.output_dir, sample + "_R1_.fastq.gz")
            fastq_r = os.path.join(self.output_dir, sample + "_R2_.fastq.gz")
            if self.option("enzyme_method") == "RAD":
                os.system("cat {} > {}".format(" ".join(self.sample_info[sample]["l"]), fastq_l))
                self.fq_list.write(sample + "\t" + fastq_l + "\n")
            else:
                os.system("cat {} > {}".format(" ".join(self.sample_info[sample]["l"]), fastq_l))
                os.system("cat {} > {}".format(" ".join(self.sample_info[sample]["r"]), fastq_r))
                self.fq_list.write(sample + "\t" + fastq_l + "\t" + fastq_r + "\n")

    def get_stacks_fastq(self):
        """
        将同一个样本同一端的fastq cat到一起
        用ngsuniform进行均一化
        RAD: 将R1端fastq作为样本的fastq
        GBS: 将R1端和R2端的fastq cat到一起作为样本的fastq
        """
        cmd_list = []
        for sample in self.sample_list:
            fastq_l = os.path.join(self.work_dir, sample + ".total.R1.fastq.gz")
            fastq_r = os.path.join(self.work_dir, sample + ".total.R2.fastq.gz")
            os.system("cat {} > {}".format(" ".join(self.sample_info[sample]["l"]), fastq_l))
            os.system("cat {} > {}".format(" ".join(self.sample_info[sample]["r"]), fastq_r))
            fastq_l_ = os.path.join(self.work_dir, sample + ".R1.fastq.gz")
            fastq_r_ = os.path.join(self.work_dir, sample + ".R2.fastq.gz")
            cmd = "{} -1 {} -2 {} -l {}".format(self.ngsuniform, fastq_l, fastq_r, self.option("uniform_length"))
            cmd += " -a {} -b {}".format(fastq_l_, fastq_r_)
            cmd_list.append(cmd)
        cmd_dict = {}
        for f in os.listdir(self.work_dir):
            if re.search("ngsuniform_cmd_|failed_ngsuniform_cmd_", f):
                os.remove(os.path.join(self.work_dir, f))
        for i in range(len(cmd_list)):
            num = i / 8
            cmd_file = os.path.join(self.work_dir, "ngsuniform_cmd_{}.list".format(str(num)))
            wrong_cmd = os.path.join(self.work_dir, "failed_ngsuniform_cmd_{}.txt".format(str(num)))
            cmd_dict[num] = [cmd_file, wrong_cmd]
            with open(cmd_file, "a") as w:
                w.write(cmd_list[i] + "\n")
        for n in cmd_dict.keys():
            cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.parafly, cmd_dict[n][0], 8, cmd_dict[n][1])
            command = self.add_command("ngsuniform_more_"+str(n), cmd_more).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("ngsuniform均一化成功")
            else:
                self.set_error("ngsuniform均一化失败，请检查", code="35501203")
        for sample in self.sample_list:
            fastq_l_ = os.path.join(self.work_dir, sample + ".R1.fastq.gz")
            fastq_r_ = os.path.join(self.work_dir, sample + ".R2.fastq.gz")
            fastq_s = os.path.join(self.output_dir, sample + ".fastq.gz")
            if os.path.exists(fastq_s):
                os.remove(fastq_s)
            if self.option("enzyme_method") == "RAD":
                os.link(fastq_l_, fastq_s)
            else:
                os.system("cat {} {} > {}".format(fastq_l_, fastq_r_, fastq_s))
            self.fq_list.write(sample + "\t" + fastq_s + "\n")

    def run(self):
        super(UniformTool, self).run()
        self.get_sample_info()
        self.fq_list = open(os.path.join(self.output_dir, "fq.list"), "w")
        if self.option("analysis_method") == "ipyrad":
            self.get_ipyrad_fastq()
        else:
            self.get_stacks_fastq()
        self.fq_list.close()
        self.end()
