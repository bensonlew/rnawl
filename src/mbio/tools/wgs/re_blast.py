# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.15

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class ReBlastAgent(Agent):
    """
    ssr比较分析:将没有比对上的scafSeq提取出来
    程序：reblast.pl
    """
    def __init__(self, parent):
        super(ReBlastAgent, self).__init__(parent)
        options = [
            {"name": "result_out", "type": "string"},  # blast_out提出的没有比对上的文件
            {"name": "scafseq_out", "type": "string"},  # scafSeq文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("result_out"):
            raise OptionError("请设置result_out文件", code="34504501")
        if not self.option("scafseq_out"):
            raise OptionError("请设置scafseq_out文件", code="34504502")
        if not os.path.exists(self.option("result_out")):
            raise OptionError("文件%s不存在，请检查",variables=(self.option("result_out")), code="34504503")
        if not os.path.exists(self.option("scafseq_out")):
            raise OptionError("文件%s不存在，请检查",variables=(self.option("scafseq_out")), code="34504504")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(ReBlastAgent, self).end()


class ReBlastTool(Tool):
    def __init__(self, config):
        super(ReBlastTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.re_blast = self.config.PACKAGE_DIR + "/wgs/reblast.pl"

    def run_re_blast(self):
        """
        运行reblast.pl
        """
        sample_name = os.path.basename(self.option("result_out")).split(".")[0]
        nomatch_seq = os.path.join(self.work_dir, sample_name + ".nomatch.seq")
        cmd = "{} {} -i {} -k {} -o {}".format(self.perl_path, self.re_blast, self.option("result_out"), self.option("scafseq_out"), nomatch_seq)
        command = self.add_command("re_blast", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("reblast.pl完成")
        else:
            self.set_error("reblast.pl失败", code="34504501")
        os.link(nomatch_seq, os.path.join(self.output_dir, sample_name + ".nomatch.seq"))

    def run(self):
        super(ReBlastTool, self).run()
        self.run_re_blast()
        self.end()
