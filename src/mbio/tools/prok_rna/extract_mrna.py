# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class ExtractMrnaAgent(Agent):
    def __init__(self, parent):
        super(ExtractMrnaAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "prok_rna.fasta"},
            {"name": "predicted_cds", "type": "infile", "format": "prok_rna.fasta"},
            {"name": "biotype", "type": "infile", "format": "prok_rna.common"},
            {"name": "out_fasta", "type": "outfile", "format": "prok_rna.fasta"},
            {"name": "rna_type", "type": "string", "default": "mRNA+sRNA"}, #mRNA or mRNA+sRNA
        ]
        self.add_option(options)
        self.step.add_steps('extract_mrna')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.extract_mrna.start()
        self.step.update()

    def step_end(self):
        self.step.extract_mrna.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("fasta").is_set:
            raise OptionError("请传入fasta序列文件")
        if not self.option("biotype").is_set:
            raise OptionError("请传入biotype文件")


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'


class ExtractMrnaTool(Tool):
    def __init__(self, config):
        super(ExtractMrnaTool, self).__init__(config)

    def extract_fasta(self):
        biotype = self.option("biotype").prop["path"]
        out_fasta = self.work_dir + "/gene.fa"
        mrna_list = {}
        with open(biotype, "r") as f:
            for line in f:
                gene_id = line.strip().split("\t")[0]
                gene_biotype = line.strip().split("\t")[1]
                if gene_biotype == "mRNA":
                    mrna_list[gene_id] = 1
        with open(out_fasta, "w") as w:
            for seq_record in SeqIO.parse(self.option('fasta').prop['path'], "fasta"):
                if self.option("rna_type") == "mRNA":
                    if seq_record.id in mrna_list:
                        w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
                    

                else:
                    if seq_record.id in mrna_list or seq_record.id.startswith("sRNA"):
                        w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
            if self.option("predicted_cds").is_set and self.option("rna_type") in ["mRNA", "mRNA+sRNA"]:
                for seq_record in SeqIO.parse(self.option("predicted_cds").prop['path'], "fasta"):
                    w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))

        self.option("out_fasta", out_fasta)

    def run(self):
        super(ExtractMrnaTool, self).run()
        self.extract_fasta()
        self.end()
