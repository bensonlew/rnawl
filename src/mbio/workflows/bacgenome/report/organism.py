# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/5/10'

from biocluster.workflow import Workflow
import os
from mbio.packages.metagenomic.common import link_file
from mbio.packages.bacgenome.common import fasta_cutoff,get_namelist


class OrganismWorkflow(Workflow):
    """
    用于完成图的Organism分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(OrganismWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "samp", "type": "string"},
            {"name": "min_len", "type": "string", "default":"1000"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sixteens = self.add_module("bacgenome.blast_gene")
        self.hgene = self.add_tool("bacgenome.hgene_blastx")

    def run(self):
        self.on_rely([self.sixteens, self.hgene], self.set_db)
        fasta_cutoff(self.option('fa').prop['path'], self.option("min_len"), self.work_dir + "/" + "seq.fa")
        self.run_sixteens()
        self.run_hgene()
        super(OrganismWorkflow, self).run()

    def run_sixteens(self):
        self.sixteens.set_options({
            "fa": self.work_dir + "/" + "seq.fa"
        })
        self.sixteens.run()

    def run_hgene(self):
        self.hgene.set_options({
            "seq_fa": self.work_dir + "/" + "seq.fa"
        })
        self.hgene.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('bac_assem.assess')  # edit api path
        self.logger.info("开始进行导表")
        if self.sixteens.option("out_table").is_set:
            sixteen_table = self.sixteens.option("out_table").prop["path"]
            api_path.add_draft_assess_16s("organism_16s", self.option("main_id"), self.option("samp"), sixteen_table)
            link_file(sixteen_table, self.output_dir + "/16s_blast.xls")
        if self.hgene.option("out_table").is_set:
            hgene_table = self.hgene.option("out_table").prop["path"]
            name_list =get_namelist(hgene_table)
            api_path.add_draft_assess_hk("organism_hk", self.option("main_id"), self.option("samp"), hgene_table, name_list)
            link_file(hgene_table, self.output_dir + "/hgene_blast.xls")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
            ["16s_blast.xls", "txt", "16s比对结果表"],
            ["hgene_blast.xls", "txt", "持家基因比对结果表"]
        ])
        super(OrganismWorkflow, self).end()
