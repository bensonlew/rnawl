# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2018.01.30

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime

starttime = datetime.datetime.now()


class GvcfTypingV2Agent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(GvcfTypingV2Agent, self).__init__(parent)
        options = [
            {"name": "vcf_list", "type": "string"},  # vcf_list文件
            {"name": "fa_file", "type": "infile", "format": "sequence.fasta"}  # 过滤后的vcf文件
        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        if not self.option("vcf_list"):
            raise OptionError("请输入vcf_list文件")
        if not self.option("fa_file"):
            raise OptionError("请输入fa_file文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(GvcfTypingV2Agent, self).end()


class GvcfTypingV2Tool(Tool):
    def __init__(self, config):
        super(GvcfTypingV2Tool, self).__init__(config)
        # 线上：10.8.0.37， 线下：10.100.201.20
        # self.set_environ(
        #     SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
        #                                                 "MajorBio_cluster_201.20.lic")   #线下
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_0.37.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        self.sentieon = "bioinfo/WGS/sentieon-genomics-201808/bin/sentieon"

    def run_gvcf_typing_v2(self):
        """
        gvcf_typing_v2

        """

        cmd = "{} driver -r {} --algo GVCFtyper".format(self.sentieon, self.option("fa_file").prop["path"])
        with open(self.option("vcf_list")) as f:
            lines = f.readlines()
            for line in lines:
                cmd += " -v {}".format(line.strip().split("\t")[1])
            cmd += " {}".format(self.output_dir + "/pop.variant.vcf")
        command = self.add_command("bam_realign", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gvcf_typing_v2运行完成")
        else:
            self.set_error("gvcf_typing_v2运行失败")

    def run(self):
        super(GvcfTypingV2Tool, self).run()
        self.run_gvcf_typing_v2()
        self.end()
