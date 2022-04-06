# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171208

"""对比对到nt数据库的xml进行统计，得到比对到物种的统计表"""
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class BlastNtStatAgent(Agent):
    def __init__(self, parent):
        super(BlastNtStatAgent, self).__init__(parent)
        options = [
            {"name": "nt_xml", "type": "infile", "format": "align.blast.blast_xml"},  # 比对到nt的xml文件
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 进行nt比对的fa序列
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("nt_xml").is_set:
            raise OptionError("必须设置比对到nt的xml")
        if not self.option("query").is_set:
            raise OptionError("必须设置进行nt比对的fasta文件")

    def set_resource(self):
        self._cpu = "1"
        self._memory = "10G"


class BlastNtStatTool(Tool):
    def __init__(self, config):
        super(BlastNtStatTool, self).__init__(config)
        self._version = 1.0
        self.perl = "program/perl-5.24.0/bin/perl"
        self.python = "miniconda2/bin/python"
        self.xml_table = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/blast_result.pl"
        self.nt_stat = self.config.PACKAGE_DIR + "/datasplit/nt_stat.py"

    def xml_to_table(self):
        """将xml文件转为table文件，没条序列只保留一条"""
        in_xml = self.option("nt_xml").prop["path"]
        out_table = os.path.basename(self.option("nt_xml").prop["path"]).split(".")[0] + ".xls"
        cmd = "{} {} -i {} -o {}".format(self.perl, self.xml_table, in_xml, out_table)
        command = self.add_command("xml_table", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("xml转table运行成功")
        else:
            self.set_error("xml转table运行失败")

    def species_stat(self):
        """统计比对到nt的table文件的物种及百分比"""
        nt_table = os.path.join(self.work_dir, os.path.basename(self.option("nt_xml").prop["path"]).split(".")[0] + ".xls")
        cmd = "{} {} -i {} -f {} -o {}".format(self.python, self.nt_stat, nt_table, self.option("query").prop["path"], "nt_species_stat.xls")
        command = self.add_command("species_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("nt物种统计运行成功")
        else:
            self.set_error("nt物种统计运行失败")
        os.link(os.path.join(self.work_dir, "nt_species_stat.xls"), os.path.join(self.output_dir, "nt_species_stat.xls"))

    def run(self):
        self.xml_to_table()
        self.species_stat()
        self.end()
