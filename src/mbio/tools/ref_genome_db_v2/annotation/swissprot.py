# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.align.blast.xml2table import *


class SwissprotAgent(Agent):
    """
    to perform Swissprot Annotation
    """
    def __init__(self, parent):
        super(SwissprotAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "swissprot_table", "type": "outfile", "format": "ref_rna_v2.blast_table"},
        ]
        self.add_option(options)
        self.step.add_steps("swissprot_anno")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.swissprot_anno.start()
        self.step.update()

    def step_end(self):
        self.step.swissprot_anno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("必须提供BLAST结果文件", code = "33703001")
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
            ["./*_vs_.*.xls", "xls", "swissprot注释结果表"]
        ])
        super(SwissprotAgent, self).end()


class SwissprotTool(Tool):
    def __init__(self, config):
        super(SwissprotTool, self).__init__(config)

    def xml2table(self):
        "xml格式转换table"
        xml_fp = self.option("blastout").prop["path"]
        table = os.path.basename(xml_fp)[:-3] + "xls"
        table_out = os.path.join(self.output_dir, table)
        xml2table(xml_fp, table_out)
        print table_out
        self.option("swissprot_table", table_out)
        self.end()

    def run(self):
        super(SwissprotTool, self).run()
        self.xml2table()
