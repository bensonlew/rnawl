# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class DiamondEnzymeAgent(Agent):
    """
    可移动元件的三种酶的比对
    """
    def __init__(self, parent):
        super(DiamondEnzymeAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "ptt", "type": "string"},  # 类似gff
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="31100201")
        return True

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(DiamondEnzymeAgent, self).end()


class DiamondEnzymeTool(Tool):
    def __init__(self, config):
        super(DiamondEnzymeTool, self).__init__(config)
        self.db_path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/database/all"
        self.db_des = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/database/all.des.xls"
        self.cmd_path = "bioinfo/align/diamond-0.9.11/diamond"
        self.python = "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/mobile_genetic_elements/enzyme_type.py"
        self.set_environ(BLASTDB=self.db_path)
        self.query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        self.outfile = self.output_dir + "/" +self.query_name + ".enzyme.xls"


    def run_blast(self):
        """
        运行diamond
        :param db_name: blastdb名称
        :return:
        """
        self.outputfile = os.path.join(self.work_dir, self.query_name)
        cmd ="{} blastp --db {} --query {} -k 1 --id 50 -o {} -e 1e-5 --threads 4".format(self.cmd_path, self.db_path, self.option("query").prop['path'], self.outputfile)
        self.logger.info("开始运行blast")
        blast_command = self.add_command("blast", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行blast完成")
        else:
            self.set_error("blast运行出错!")

    def run_get_type(self):
        cmd ="{} {} --in {} --db {} --gff {} --o {}".format(self.python, self.python_script, self.outputfile, self.db_des, self.option("ptt"), self.outfile)
        self.logger.info("开始运行get_type")
        blast_command = self.add_command("get_type", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行get_type完成")
        else:
            self.set_error("get_type运行出错!")

    def run(self):
        """
        运行
        :return:
        """
        super(DiamondEnzymeTool, self).run()
        self.run_blast()
        self.run_get_type()
        self.end()