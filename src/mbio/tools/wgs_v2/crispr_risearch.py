# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# modified 20190320

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class CrisprRisearchAgent(Agent):
    """
    接口“变异位点比较分析Tool, 导表在tool中完成”
    """

    def __init__(self, parent):
        super(CrisprRisearchAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "string"},
            {"name": "sgrna_fa", "type": "infile", "format": "wgs_v2.fasta"},
        ]
        self.add_option(options)
        self.step.add_steps('crispr_riserch')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.crispr_riserch.start()
        self.step.update()

    def step_end(self):
        self.step.crispr_riserch.finish()
        self.step.update()

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("必须输入ref.fa")
        if not self.option("sgrna_fa"):
            raise OptionError("必须输入sgrna.fa文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(CrisprRisearchAgent, self).end()


class CrisprRisearchTool(Tool):

    def __init__(self, config):
        super(CrisprRisearchTool, self).__init__(config)
        self.riserch = 'bioinfo/wgs_v2/risearch2.x'
        self.perl = "program/perl/perls/perl-5.24.0/bin/perl"
        self.sgrna_region = self.config.PACKAGE_DIR + "/wgs_v2/sgrna.region.pl"

    def run_riserch2(self):

        cmd = "{} -c {} -o {}".format(self.riserch, os.path.join(self.config.SOFTWARE_DIR + "/database/dna_geneome",
                                                                 self.option("ref_fa")), self.output_dir + "/ref.suf")
        # cmd = "{} -c {} -o {}".format(self.riserch, self.option("ref_fa"), self.output_dir + "/ref.suf")
        command = self.add_command("riserch_step1", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ref.suf生成完成")
        else:
            self.set_error("ref.suf生成失败")
        cmd2 = "{} -q {} -i {} -s 1:20 -m 4:1 -e 10000 -l 0 --noGUseed -p3 -t 8".format(self.riserch,
                                                                                        self.option("sgrna_fa").
                                                                                        prop["path"],
                                                                                        self.output_dir + "/ref.suf")
        command1 = self.add_command("riserch_step2", cmd2).run()
        self.wait()
        if command1.return_code == 0:
            self.logger.info("riserch运行完成")
        else:
            self.set_error("riserch运行失败")

    def run_sgrna_region(self):
        """"
        用于生成下一步要使用的region文件"""
        cmd = "{} {} -i {} -o {} -sgr {}".format(self.perl, self.sgrna_region, self.work_dir +
                                                 "/risearch_sgrna.out.gz", self.output_dir + "/sgRNA.region",
                                                 self.option("sgrna_fa").prop["path"])
        command = self.add_command("riserch_step3", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("sgRNA区域生成完成")
        else:
            self.set_error("sgRNA区域生成失败")

    def run(self):
        super(CrisprRisearchTool, self).run()
        self.run_riserch2()
        self.run_sgrna_region()
        self.end()
