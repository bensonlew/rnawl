# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.03.12

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class GfatofastaAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(GfatofastaAgent, self).__init__(parent)
        options = [
            {"name": "gfa_file", "type": "infile", "format": "tool_lab.no_empty"},
            {"name": "sample_name", "type": "string"},  # 样本名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gfa_file").is_set:
            raise OptionError("请设置gfa_file")
        if not self.option("sample_name"):
            raise OptionError("请设置sample_name")

    def set_resource(self):
        self._cpu = 6
        self._memory = "20G"

    def end(self):
        super(GfatofastaAgent, self).end()


class GfatofastaTool(Tool):
    def __init__(self, config):
        super(GfatofastaTool, self).__init__(config)
        self.sh = 'awk \'/^S/{print ">"$2"\\n"$3}\' '

    def run_awk(self):
        """

        """

        os.system("chmod 777 {}".format(self.option("gfa_file").prop["path"]))
        self.logger.info(self.sh + ' {} > {}'.format(self.option("gfa_file").prop["path"], self.output_dir + "/" + self.option("sample_name") +".fasta"))
        return_code = os.system(self.sh + ' {} > {}'.format(self.option("gfa_file").prop["path"], self.output_dir + "/" + self.option("sample_name") +".fasta"))
        #command = self.add_command("run_awk", cmd).run()
        self.wait()
        if return_code == 0:
            self.logger.info("run_awk运行成功")
        else:
            self.set_error("run_awk运行失败")

    def run(self):
        super(GfatofastaTool, self).run()
        self.run_awk()
        self.end()
