# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __modify__ = '2020/2/10'

from biocluster.workflow import Workflow
import os
from mbio.packages.metagenomic.common import link_file,get_namelist


class CompleteOrganismWorkflow(Workflow):
    """

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CompleteOrganismWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "genome_id", "type": "string"},
            {"name": "samp", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sixteens = self.add_module("bacgenome.blast_gene")
        self.hgene = self.add_tool("bacgenome.hgene_blastx")

    def run(self):
        self.on_rely([self.sixteens, self.hgene], self.set_db)
        self.run_sixteens()
        self.run_hgene()
        super(CompleteOrganismWorkflow, self).run()

    def run_sixteens(self):
        self.sixteens.set_options({
            "fa": self.option("fa")
        })
        self.sixteens.run()

    def run_hgene(self):
        self.hgene.set_options({
            "seq_fa": self.option("fa")
        })
        self.hgene.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('bac_assem.blast')  # edit api path
        self.logger.info("开始进行导表")
        if self.sixteens.option("out_table").is_set:
            sixteen_table = self.sixteens.option("out_table").prop["path"]
            api_path.add_complete_organism_16s("complete_organism_16s", self.option("main_id"), self.option("genome_id"), self.option("samp"),
                                      sixteen_table)
            link_file(sixteen_table, self.output_dir + "/16s_blast.xls")
        if self.hgene.option("out_table").is_set:
            hgene_table = self.hgene.option("out_table").prop["path"]
            name_list = get_namelist(hgene_table)
            api_path.add_complete_organism_hk("complete_organism_hk", self.option("main_id"), self.option("genome_id"), self.option("samp"), hgene_table, name_list)
            link_file(hgene_table, self.output_dir + "/hgene_blast.xls")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
            ["16s_blast.xls", "txt", "16s比对结果表"],
            ["hgene_blast.xls", "txt", "持家基因比对结果表"]
        ])
        super(CompleteOrganismWorkflow, self).end()