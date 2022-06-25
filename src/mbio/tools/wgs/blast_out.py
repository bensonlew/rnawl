# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.08

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class BlastOutAgent(Agent):
    """
    ssr比较分析处理blastn后的结果
    程序：blastout.pl
    """
    def __init__(self, parent):
        super(BlastOutAgent, self).__init__(parent)
        options = [
            {"name": "blast_out", "type": "string"},  # blastn输出结果文件
            {"name": "scafseq_out", "type": "string"},  # scafSeq文件
            {"name": "ssr_result", "type": "string"},  # ssr结果文件result.xls
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("blast_out"):
            raise OptionError("请设置blast_out文件", code="34500701")
        if not self.option("scafseq_out"):
            raise OptionError("请设置scafseq_out文件", code="34500702")
        if not self.option("ssr_result"):
            raise OptionError("请设置ssr_result文件", code="34500703")
        if not os.path.exists(self.option("blast_out")):
            raise OptionError("文件%s不存在，请检查", variables=(self.option("blast_out")), code="34500704")
        if not os.path.exists(self.option("scafseq_out")):
            raise OptionError("文件%s不存在，请检查", variables=(self.option("scafseq_out")), code="34500705")
        if not os.path.exists(self.option("ssr_result")):
            raise OptionError("文件%s不存在，请检查", variables=(self.option("ssr_result")), code="34500706")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(BlastOutAgent, self).end()


class BlastOutTool(Tool):
    def __init__(self, config):
        super(BlastOutTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.blast_out = self.config.PACKAGE_DIR + "/wgs/blastout.pl"

    def run_blast_out(self):
        """
        运行blastout.pl
        """
        sample_name = os.path.basename(self.option("blast_out")).split(".")[0]
        scaf_seq = os.path.join(self.work_dir, os.path.basename(self.option("scafseq_out")))
        os.link(self.option("scafseq_out"), scaf_seq)
        scafseq_path = scaf_seq.split(".")[0]
        cmd = "{} {} -i {} -o {} -k {}".format(self.perl_path, self.blast_out, self.option("blast_out"), scafseq_path,
                                               self.option("ssr_result"))
        command = self.add_command("blast_out", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("blastout.pl完成")
        else:
            self.set_error("blastout.pl失败", code="34500701")
        os.link(os.path.join(self.work_dir, sample_name + ".result-1.out"), os.path.join(self.output_dir, sample_name + ".result-1.out"))
        os.link(os.path.join(self.work_dir, sample_name + ".result-2.out"), os.path.join(self.output_dir, sample_name + ".result-2.out"))

    def run(self):
        super(BlastOutTool, self).run()
        self.run_blast_out()
        self.end()
