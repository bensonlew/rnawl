# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class ExtractPttAgent(Agent):
    def __init__(self, parent):
        super(ExtractPttAgent, self).__init__(parent)
        options = [
            {"name": "ptt", "type": "infile", "format": "prok_rna.common"},
            {"name": "biotype", "type": "infile", "format": "prok_rna.common"},
            {"name": "out_ptt", "type": "outfile", "format": "prok_rna.common"},
            {"name": "rna_type", "type": "string", "default": "mRNA+sRNA"}, #mRNA or mRNA+sRNA
        ]
        self.add_option(options)
        self.step.add_steps('extract_ptt')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.extract_ptt.start()
        self.step.update()

    def step_end(self):
        self.step.extract_ptt.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("ptt").is_set:
            raise OptionError("请传入ptt.bed文件")
        if not self.option("biotype").is_set:
            raise OptionError("请传入biotype文件")


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '3G'


class ExtractPttTool(Tool):
    def __init__(self, config):
        super(ExtractPttTool, self).__init__(config)

    def extract_ptt(self):
        ptt = self.option("ptt").prop["path"]
        biotype = self.option("biotype").prop["path"]
        out_ptt = self.work_dir + "/ptt.bed"
        mrna_list = {}
        with open(biotype, "r") as f:
            for line in f:
                gene_id = line.strip().split("\t")[0]
                gene_biotype = line.strip().split("\t")[1]
                if gene_biotype == "mRNA":
                    mrna_list[gene_id] = 1
        with open (ptt, "r") as f, open(out_ptt, "w") as w:
            for line in f:
                items = line.strip().split("\t")
                if items[6] in mrna_list:
                    w.write(line)
        self.option("out_ptt", out_ptt)

    def run(self):
        super(ExtractPttTool, self).run()
        self.extract_ptt()
        self.end()
