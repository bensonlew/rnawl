# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.03.12

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class InsertResultAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(InsertResultAgent, self).__init__(parent)
        options = [
            {"name": "result_csv", "type": "infile", "format": "wgs_v2.bcf"},  # 279_1.result.csv
            {"name": "sample", "type": "string"}  # 样本名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("result_csv").is_set:
            raise OptionError("请设置result_csv")
        if not self.option("sample"):
            raise OptionError("请设置sample")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(InsertResultAgent, self).end()


class InsertResultTool(Tool):
    def __init__(self, config):
        super(InsertResultTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.insertresult = "/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgs_v2/insertresult.pl"

    def run_insertresult(self):
        """

        """
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.insertresult, self.option("result_csv").prop["path"],
                                         os.path.join(self.output_dir, (self.option("sample") + ".insert.xls")))
        command = self.add_command("insertresult", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("insertresult运行成功")
        else:
            self.set_error("insertresult运行失败")

    def run(self):
        super(InsertResultTool, self).run()
        self.run_insertresult()
        self.end()
