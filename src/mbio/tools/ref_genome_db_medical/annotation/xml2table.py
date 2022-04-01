# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.align.blast.xml2table import *


class Xml2tableAgent(Agent):
    """
    将注释的xml文件转为table文件
    """
    def __init__(self, parent):
        super(Xml2tableAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast_table", "type": "outfile", "format": "ref_rna_v2.blast_table"},
        ]
        self.add_option(options)
        self.step.add_steps("xml2table")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.xml2table.start()
        self.step.update()

    def step_end(self):
        self.step.xml2table.finish()
        self.step.update()

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("必须提供BLAST结果文件", code = "33703101")
        else:
            pass

    def set_resource(self):
        self._cpu = 10
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["./*_vs_.*.xls", "xls", "blast结果表"]
        ])
        super(Xml2tableAgent, self).end()


class Xml2tableTool(Tool):
    def __init__(self, config):
        super(Xml2tableTool, self).__init__(config)

    def xml2table(self):
        "xml格式转换table"
        xml_fp = self.option("blastout").prop["path"]
        table = os.path.basename(xml_fp)[:-3] + "xls"
        table_out = os.path.join(self.output_dir, table)
        xml2table(xml_fp, table_out)
        self.option("blast_table", table_out)
        self.end()

    def run(self):
        super(Xml2tableTool, self).run()
        self.xml2table()
