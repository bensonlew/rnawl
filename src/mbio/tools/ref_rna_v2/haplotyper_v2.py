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
            {"name": "bam_file", "type": "infile", "format": "ref_rna_v2.common"},  # 过滤后的vcf文件
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},  # 过滤后的vcf文件
            {"name": "name", "type": "string"}  # 生成文件名字，测试文件中名字为DE1_10.g.vcf，其中.g.vcf为固定。

        ]
        self.add_option(options)
        self._memory_increase_step = 80

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("请设置bam路径", code="33709503")
        if not self.option("fa_file"):
            raise OptionError("请设置ref.fa路径", code="33709504")

    def set_resource(self):
        self._cpu = 10
        self._memory = "100G"

    def end(self):
        super(HaplotyperV2Agent, self).end()


class HaplotyperV2Tool(Tool):
    def __init__(self, config):
        super(HaplotyperV2Tool, self).__init__(config)
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_201.20.lic")
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
        command = self.add_command("bam_realign", cmd, ignore_error=True).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("haplotyper_v2运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("haplotyper_v2运行失败", code="33709502")

    def run(self):
        super(HaplotyperV2Tool, self).run()
        self.run_haplotyper_v2()
        self.end()
