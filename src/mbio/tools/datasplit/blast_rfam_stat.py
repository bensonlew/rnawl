# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171215

"""对比对到Rfam数据库的xml文件进行统计，得到比对到Rfam的分类统计表"""
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class BlastRfamStatAgent(Agent):
    def __init__(self, parent):
        super(BlastRfamStatAgent, self).__init__(parent)
        options = [
            {"name": "rfam_xml", "type": "infile", "format": "align.blast.blast_xml"},  # 比对到Rfam的xml文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("rfam_xml").is_set:
            raise OptionError("必须设置比对到Rfam的xml文件")

    def set_resource(self):
        self._cpu = "1"
        self._memory = "5G"

class BlastRfamStatTool(Tool):
    def __init__(self, config):
        super(BlastRfamStatTool, self).__init__(config)
        self._version = 1.0
        self.python = "miniconda2/bin/python"
        self.rfam_stat = self.config.PACKAGE_DIR + "/datasplit/blast_rfam_stat.py"
        # self.rfam_seed = self.config.SOFTWARE_DIR + "/database/Rfam/Rfam.seed"
        self.rfam_seed = self.config.SOFTWARE_DIR + "/database/align/ncbi/db/rfam_v14.6/Rfam.seed"

    def rna_taxon_stat(self):
        """统计比对到Rfam的xml文件对应Rfam.seed的RNA分类"""
        rfam_xml = self.option("rfam_xml").prop["path"]
        summary_out = self.output_dir + "/rfam_summary.xls"
        cmd = "{} {} -i {} -db {} -o {}".format(self.python, self.rfam_stat, rfam_xml, self.rfam_seed, summary_out)
        command = self.add_command("rfam_summary", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("Rfam统计运行完成")
        else:
            self.set_error("Rfam统计运行失败")

    def run(self):
        self.rna_taxon_stat()
        self.end()
