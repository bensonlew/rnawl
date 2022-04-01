# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.02.20

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class DellyCallAgent(Agent):
    """
    工具：delly-call
    """
    def __init__(self, parent):
        super(DellyCallAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "string"},  # ref.fa
            {"name": "bam_list", "type": "string"},  # 最终的bam文件路径
            {"name": "bcf_file", "type": "outfile", "format": "wgs_v2.bcf"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa")
        if not self.option("bam_list"):
            raise OptionError("请设置bam_list")

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(DellyCallAgent, self).end()


class DellyCallTool(Tool):
    def __init__(self, config):
        super(DellyCallTool, self).__init__(config)
        self.delly = "bioinfo/wgs_v2/delly"
        # self.ref_fa = os.path.join(self.config.SOFTWARE_DIR, ("database/dna_geneome/" + self.option("ref_fa")))

    def run_delly_call(self):
        """
        delly call -t ALL -g ref.fa -o xxx.bcf sample1.bam sample2.bam
        """
        output_path = os.path.join(self.output_dir, "pop.nosort.sv.bcf")
        cmd = "{} call -t All -g {} -o {}".format(self.delly, self.option("ref_fa"), output_path)
        with open(self.option("bam_list"), "r")as fr:
            lines = fr.readlines()
            for line in lines:
                tmp = line.strip().split("\t")
                cmd += " {}".format(tmp[1])
        command = self.add_command("delly_call", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("delly_call运行成功")
        else:
            self.set_error("delly_call运行失败")
        self.option("bcf_file", output_path)

    def run(self):
        super(DellyCallTool, self).run()
        self.run_delly_call()
        self.end()
