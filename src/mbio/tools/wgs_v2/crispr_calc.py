# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# modified 20190320

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class CrisprCalcAgent(Agent):
    """
    接口“变异位点比较分析Tool, 导表在tool中完成”
    """

    def __init__(self, parent):
        super(CrisprCalcAgent, self).__init__(parent)
        options = [
            {"name": "wt_dir", "type": "infile",
             "format": "wgs_v2.wt_dir"},
            {"name": "sgrna_region", "type": "infile",
             "format": "wgs_v2.sgrna_region"},
            {"name": "mut_dir", "type": "infile", "format": "wgs_v2.wt_dir"},
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('crispr_calc')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.crispr_calc.start()
        self.step.update()

    def step_end(self):
        self.step.crispr_calc.finish()
        self.step.update()

    def check_options(self):
        if not self.option("wt_dir"):
            raise OptionError("必须输入野生型文件")
        if not self.option("mut_dir"):
            raise OptionError("必须输入突变型文件")
        if not self.option("sgrna_region"):
            raise OptionError("必须输入region文件")
        if not self.option("name"):
            raise OptionError("必须输入生成文件名称")

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(CrisprCalcAgent, self).end()


class CrisprCalcTool(Tool):

    def __init__(self, config):
        super(CrisprCalcTool, self).__init__(config)
        self.perl = "miniconda2/bin/perl"
        self.crispr_stat = self.config.PACKAGE_DIR + "/wgs_v2/crispr.pl"

    def run_crispr_stat(self):
        cmd = "{} {} -wt {} -mut {} -region {} -output {}".format(
            self.perl,
            self.crispr_stat,
            self.option("wt_dir").prop["path"],
            self.option("mut_dir").prop["path"],
            self.option("sgrna_region").prop["path"],
            self.output_dir + "/" + self.option("name"))
        command = self.add_command("crispr_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("crispr_stat 计算完成")
        else:
            self.set_error("crispr_stat 计算失败")


    # def run_calc_mut(self):
    #     """
    #     用于统计变异类型
    #     :return:
    #     """
    #     R_num = ''
    #     chr = self.option("region").strip().split(":")[0]
    #     span = range(int(self.option("region").strip().split(":")[1].strip().split("-")[0]),
    #                  int(self.option("region").strip().split(":")[1].strip().split("-")[1])
    #                  + 1)
    #     with open(self.option("sgrna_region").prop["path"]) as f:
    #         # with open(self.work_dir + "/CrisprRisearch/output/sgRNA.region")
    #         # as f:
    #         lines = f.readlines()
    #         for line in lines:
    #             chr_f = line.strip().split("\t")[0]
    #             if chr_f == chr:
    #                 if int(
    #                     line.strip().split("\t")[1]) in span and int(
    #                         line.strip().split("\t")[2]) in span:
    #                     R_num = line.strip().split("\t")[3]
    #     print "#################################"
    #     print R_num
    #     print "#################################"
    #     path = self.option("mut_dir").prop["path"]
    #     samples = os.listdir(path)
    #     sum_insertions = 0
    #     sum_deletions = 0
    #     sum_substitutions = 0
    #     with open(self.output_dir + "/calculations", "w") as w:
    #         w.write("insertions\tdeletions\tsubstitutions\n")
    #         for sample in samples:
    #             if os.path.exists(
    #                 path +
    #                 "/" +
    #                 sample +
    #                 "/CRISPResso_on_" +
    #                 R_num +
    #                 "/"
    #                 "Quantification_of_editing_freque"
    #                     "ncy.txt"):
    #                 with open(path + "/" + sample + "/CRISPResso_on_" + R_num + "/Quantification_of_editing_"
    #                                                                             "frequency.txt") as f:
    #                     lines = f.readlines()
    #                     insertions = 0
    #                     deletions = 0
    #                     substitutions = 0
    #                     for i in range(2, 5):
    #                         insertions += int(
    #                             ",".join(
    #                                 re.findall(
    #                                     ".*reads \((.*) reads with insertions.*", lines[i])))
    #                         deletions += int(
    #                             ",".join(
    #                                 re.findall(
    #                                     ".*insertions, (.*) reads with deletions.*",
    #                                     lines[i])))
    #                         substitutions += int(
    #                             ",".join(
    #                                 re.findall(
    #                                     ".*deletions, (.*) reads with substitutions.*",
    #                                     lines[i])))
    #                     print insertions
    #                     print deletions
    #                     print substitutions
    #                     print "++++++++++++++++++++++++++++++++++++++++++++++"
    #                     sum_insertions += insertions
    #                     sum_deletions += deletions
    #                     sum_substitutions += substitutions
    #             else:
    #                 w.write("## None\n")
    #         w.write(
    #             str(sum_insertions) +
    #             "\t" +
    #             str(sum_deletions) +
    #             "\t" +
    #             str(sum_substitutions) +
    #             "\n")
    #         print sum_deletions
    #         print sum_deletions
    #         print sum_substitutions
    #         print "*88888888888888888888888888888888"
    #         self.logger.info("88888888888888888888888888888888")

    def run(self):
        super(CrisprCalcTool, self).run()
        self.run_crispr_stat()
        self.end()
