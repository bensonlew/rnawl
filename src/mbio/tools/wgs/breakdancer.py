# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.12

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class BreakdancerAgent(Agent):
    """
    软件:breakdancer
    call sv
    """
    def __init__(self, parent):
        super(BreakdancerAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.sam"},  # bam文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file").is_set:
            raise OptionError("请设置sam文件", code="34500801")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(BreakdancerAgent, self).end()


class BreakdancerTool(Tool):
    def __init__(self, config):
        super(BreakdancerTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/align/samtools-1.7")
        self.perl_path = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin/perl"
        self.bam2cfg_sh = "bioinfo/WGS/breakdancer_bam2cfg.sh"
        self.bam2cfg_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/breakdancer-1.3.6/perl/bam2cfg.pl"
        self.breakdancer_max_sh = "bioinfo/WGS/breakdancer_max.sh"
        self.breakdancer_max = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/breakdancer-1.3.6/build/bin/breakdancer-max"

    def run_bam2cfg(self):
        """
        bam2cfg.pl
        """
        self.sample_name = os.path.basename(self.option("bam_file").prop["path"]).split(".")[0]
        cfg_file = os.path.join(self.work_dir, self.sample_name + ".cfg")
        cmd = "{} {} {} 5 {} {}".format(self.bam2cfg_sh, self.perl_path, self.bam2cfg_path, self.option("bam_file").prop["path"], cfg_file)
        command = self.add_command("sv_bam2cfg", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam2cfg.pl运行完成")
        else:
            self.set_error("bam2cfg.pl运行失败", code="34500801")
        if os.path.getsize(cfg_file) == 0:
            self.logger.info("变异系数大于截断值，重新调整截断值进行分析")
            cfg_cmd = os.path.join(self.work_dir, "sv_bam2cfg.o")
            with open(cfg_cmd, "r") as f:
                lines = f.readlines()
                line = lines[-2].strip().split(" ")
                cutoff = float(line[3])
                cutoff = int(cutoff) + 1
            self.logger.info("截断值：{}".format(cutoff))
            cmd = "{} {} {} {} {} {}".format(self.bam2cfg_sh, self.perl_path, self.bam2cfg_path, cutoff, self.option("bam_file").prop["path"], cfg_file)
            command = self.add_command("sv_bam2cfg1", cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("bam2cfg.pl运行完成")
            else:
                self.set_error("bam2cfg.pl运行失败", code="34500802")
            if os.path.getsize(cfg_file) == 0:
                self.set_error("没有生成cfg文件，截断值可能需要调整，请检查", code="34500803")

    def run_breakdancer_max(self):
        """
        breakdancer-max
        """
        cmd = "{} {} 50 10 {} {} {}".format(self.breakdancer_max_sh, self.breakdancer_max, self.work_dir + "/" + self.sample_name + ".ctx",
                                            os.path.join(self.work_dir, self.sample_name + ".cfg"), os.path.join(self.work_dir, self.sample_name + ".sv"))
        command = self.add_command("breakdancer_max", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("breakdancer_max运行完成")
        else:
            self.set_error("breakdancer_max运行失败", code="34500804")
        if os.path.exists(os.path.join(self.output_dir, self.sample_name + ".sv.xls")):
            os.remove(os.path.join(self.output_dir, self.sample_name + ".sv.xls"))
        os.link(os.path.join(self.work_dir, self.sample_name + ".sv"), os.path.join(self.output_dir, self.sample_name + ".sv.xls"))

    def run(self):
        super(BreakdancerTool, self).run()
        self.run_bam2cfg()
        self.run_breakdancer_max()
        self.end()
