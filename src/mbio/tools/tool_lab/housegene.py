#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'zzg'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os, shutil


class HousegeneAgent(Agent):
    """
    持家基因
    version 1.0
    author: zzg
    last_modify: 2021.2.20
    """

    def __init__(self, parent):
        super(HousegeneAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('fasta').is_set:
            raise OptionError("请传入序列文件！", code="")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 8
        self._memory = '10G'

    def end(self):
        super(HousegeneAgent, self).end()


class HousegeneTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(HousegeneTool, self).__init__(config)
        self.predict_seq = self.add_module("bac_comp_genome.gene_predict_seq")
        self.core_gene = self.add_tool('bac_comp_genome.get_core_gene')

    def run_gene_predict(self):
        opts = ({
            "genome": self.option("fasta"),
            "genome_name": self.option("sample_name"),
            "gene_prefix": "gene",
            "genome_type": "draft"
        })
        self.predict_seq.set_options(opts)
        self.predict_seq.run()

    def run_coregene(self):
        opts = {
            "seq_faa": self.predict_seq.output_dir + "/" + self.option("sample_name") + "_CDS.faa",
            "seq_gff": self.predict_seq.output_dir + "/" + self.option("sample_name") + "_CDS.gff",
            "sample_name": self.option("sample_name"),
            "method": "genome",
        }
        self.core_gene.set_options(opts)
        self.core_gene.run()

    def set_output(self):
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + "/" + self.option("sample_name")):
            shutil.rmtree(self.output_dir + "/" + self.option("sample_name"))
        os.mkdir(self.output_dir + "/" + self.option("sample_name"))
        os.link(self.core_gene.output_dir + '/' + self.option("sample_name") + '.cor_gene.fa',
                self.output_dir + "/" + self.option("sample_name") + '/' + self.option("sample_name") + '_all.hgene.fa')
        raw_hgene_file = self.core_gene.output_dir + '/' + self.option("sample_name") + '.result.xls'
        hgenefa_file = self.output_dir + "/" + self.option("sample_name") + '/' + self.option(
            "sample_name") + '_hgene.fa'
        blast_file = self.output_dir + "/" + self.option("sample_name") + '/' + self.option(
            "sample_name") + '_hgene_blast.xls'
        stat_file = self.output_dir + "/" + 'stat_file.txt'
        with open(raw_hgene_file, "r") as r, open(hgenefa_file, "w") as h, open(blast_file, "w") as b, open(stat_file,
                                                                                                            "w") as s:
            data = r.readlines()
            s.write(self.option("sample_name") + "\t" + str(len(data) - 1))
            b.write(
                "Genome" + "\t" + "Gene ID" + "\t" + "Location" + "\t" + "Start" + "\t" + "End" + "\t" + "Name" + "\t" + "Indentity" + "\t" + "Coverage" + "\n")
            for i in data[1:]:
                b.write("\t".join(i.split("\t")[0:8]) + "\n")
                h.write(">" + i.split("\t")[0] + "|" + i.split("\t")[1] + "|" + i.split("\t")[2] + ":" + i.split("\t")[
                    3] + "-" + i.split("\t")[4] + "|" + i.split("\t")[5] + "\n" + i.split("\t")[8])
        if os.path.getsize(hgenefa_file) == 0:
            shutil.rmtree(self.output_dir + "/" + self.option("sample_name"))
        self.end()

    def run(self):
        """
        运行
        """
        super(HousegeneTool, self).run()
        self.run_gene_predict()
        self.run_coregene()
        self.set_output()
        self.end()