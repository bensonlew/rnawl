# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.04.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class BwaIndexAgent(Agent):
    """
    软件:bwa index
    """
    def __init__(self, parent):
        super(BwaIndexAgent, self).__init__(parent)
        options = [
            {"name": "insert_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("insert_fa"):
            raise OptionError("必须输入insert.fa文件", code="34501001")
        if not self.option("ref_fa"):
            raise OptionError("必须输入ref.fa文件", code="34501002")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(BwaIndexAgent, self).end()


class BwaIndexTool(Tool):
    def __init__(self, config):
        super(BwaIndexTool, self).__init__(config)
        self.bwa_path = 'bioinfo/align/bwa-0.7.15/bwa'

    def run_cat(self):
        """
        cat
        """
        r2 = os.path.join(self.work_dir, 'insert_new.fa')
        with open(self.option("insert_fa").prop["path"], 'r')as f1:
            with open(r2, 'w') as fp:
                lines = f1.readlines()
                fa_num = 0
                for line in lines:
                    match_fa = re.match("^>", line)
                    if match_fa:
                        line = '>Chrinsert' + '\n'
                        fa_num += 1
                    fp.write(line)
                if fa_num < 1:
                    self.set_error("序列少于1条", code="34501001")
                elif fa_num > 1:
                    self.set_error("序列大于1条", code="34501002")
        code = os.system('cat {} {} > {}/pop.fa'.format(r2, self.option("ref_fa").prop["path"], self.output_dir))
        if code == 0:
            self.logger.info("cat成功！")
        else:
            self.set_error("cat失败！", code="34501003")

    def run_bwa(self):
        """
        bwa index
        """
        cmd = "{} index {}".format(self.bwa_path, self.output_dir + "/pop.fa")
        command = self.add_command("fa_bwa", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bwa运行完成")
        else:
            self.set_error("bwa运行失败", code="34501004")

    def run(self):
        super(BwaIndexTool, self).run()
        self.run_cat()
        self.run_bwa()
        self.end()
