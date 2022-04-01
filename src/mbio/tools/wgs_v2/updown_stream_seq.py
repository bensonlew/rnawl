# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.03.12

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class UpdownStreamSeqAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(UpdownStreamSeqAgent, self).__init__(parent)
        options = [
            {"name": "insert_xls", "type": "infile", "format": "wgs_v2.bcf"},  # 279_1.result.csv
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # ref.fa
            {"name": "sample", "type": "string"}  # 样本名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("insert_xls").is_set:
            raise OptionError("请设置insert_xls")
        if not self.option("ref_fa").is_set:
            raise OptionError("请设置ref_fa")
        if not self.option("sample"):
            raise OptionError("请设置sample")

    def set_resource(self):
        self._cpu = 4
        self._memory = "50G"

    def end(self):
        super(UpdownStreamSeqAgent, self).end()


class UpdownStreamSeqTool(Tool):
    def __init__(self, config):
        super(UpdownStreamSeqTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.updownstreamseq = "/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgs_v2/updownstreamseq.pl"

    def run_updownstreamseq(self):
        """

        """
        outfile = os.path.join(self.output_dir, self.option("sample"))
        cmd = "{} {} -i {} -fa {} -o {}".format(self.perl_path, self.updownstreamseq,
                                                self.option("insert_xls").prop["path"],
                                                self.option("ref_fa").prop["path"], outfile)
        command = self.add_command("updownstreamseq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("updownstreamseq运行成功")
        else:
            self.set_error("updownstreamseq运行失败")

    def run(self):
        super(UpdownStreamSeqTool, self).run()
        self.run_updownstreamseq()
        self.end()
