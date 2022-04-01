# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2018.12.13

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class HaplotyperV2Agent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(HaplotyperV2Agent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.bam"},  # 过滤后的vcf文件
            {"name": "fa_file", "type": "infile", "format": "sequence.fasta"},  # 过滤后的vcf文件
            {"name": "name", "type": "string"}  # 生成文件名字，测试文件中名字为DE1_10.g.vcf，其中.g.vcf为固定。

        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("请设置bam路径")
        if not self.option("fa_file"):
            raise OptionError("请设置ref.fa路径")

    def set_resource(self):
        self._cpu = 2
        self._memory = "30G"

    def end(self):
        super(HaplotyperV2Agent, self).end()


class HaplotyperV2Tool(Tool):
    def __init__(self, config):
        super(HaplotyperV2Tool, self).__init__(config)
        # 线上：10.8.0.37， 线下：10.100.201.20
        # self.set_environ(
        #     SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
        #                                                 "MajorBio_cluster_201.20.lic")    # 线下
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_0.37.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        self.sentieon = "bioinfo/WGS/sentieon-genomics-201808/bin/sentieon"

    def run_haplotyper_v2(self):
        """
        PopLDdecay

        """
        cmd = "{} driver -t 8 -i {} -r {} --algo Haplotyper --emit_mode GVCF {}"\
            .format(self.sentieon, self.option("bam_file").prop["path"], self.option("fa_file").prop["path"],
                    self.output_dir + "/" + self.option("name"))
        command = self.add_command("bam_realign", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("haplotyper_v2运行完成")
        else:
            self.set_error("haplotyper_v2运行失败")

    def run(self):
        super(HaplotyperV2Tool, self).run()
        self.run_haplotyper_v2()
        self.end()
