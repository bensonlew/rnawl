# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2018.09.12

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime

starttime = datetime.datetime.now()


class LdDecayAgent(Agent):
    """
    连锁不平衡接口Tool
    """

    def __init__(self, parent):
        super(LdDecayAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf"},  # 过滤后的vcf文件
            {"name": "maf", "type": "float"},
            {"name": "miss", "type": "float", "default": 0.3},
            {"name": "group_name", "type": "string"},  # 用于标记group的名字，生成的结果命名用
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("vcf_file"):
            raise OptionError("请设置vcf.file路径", code="34901001")

    def set_resource(self):
        self._cpu = 2
        self._memory = "13G"

    def end(self):
        super(LdDecayAgent, self).end()


class LdDecayTool(Tool):
    def __init__(self, config):
        super(LdDecayTool, self).__init__(config)
        self.R_path = 'program/R-3.3.1/bin/Rscript'
        self.ld_decay = self.config.PACKAGE_DIR + "/dna_evolution/ld-decay.R"
        self.PopLDdecay = "bioinfo/dna_evolution/PopLDdecay"

    def run_pop_ld_decay(self):
        """
        PopLDdecay

        """
        with open(self.option("vcf_file").prop["path"]) as f, open(os.path.join(self.work_dir, "new_vcf"), "w") as w:
            # 生成新的vcf文件，没有前面以双##开头的注释。
            for line in f:
                if line.startswith("##"):
                    pass
                else:
                    w.write(line)
        cmd1 = "{} --InVCF {} --OutStat {} --MAF {} --Miss {}".format(
            self.PopLDdecay,
            os.path.join(self.work_dir, "new_vcf"),
            os.path.join(self.output_dir, self.option("group_name")),
            self.option("maf"),
            1 - self.option("miss"))

        command = self.add_command("pop_ld_decay", cmd1).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("PopLDdecay运行完成")
        else:
            self.set_error("PopLDdecay运行失败", code="34901001")
        # cmd2 = "{} {} --infile {} --outfile {}".format(                        这个好像是用于画图的，可能不需要
        #     self.R_path, self.ld_decay, os.path.join(
        #         self.output_dir, "pop1.stat.gz"), os.path.join(
        #         self.output_dir, "pop1.ld"))  # 这里的pop1是需要调整的。
        # command = self.add_command("depth_stat_windows", cmd2).run()
        # self.wait()
        # if command.return_code == 0:
        #     self.logger.info("R脚本ld_decay运行完成")
        # else:
        #     self.set_error("R脚本ld_decay运行失败")

    def run(self):
        super(LdDecayTool, self).run()
        self.run_pop_ld_decay()
        self.end()
