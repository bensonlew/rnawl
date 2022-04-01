# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.09.24

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class AnnotationModule(Module):
    """
    单个细菌基因组数据库注释模块；分析内容：cog,kegg,cazy,card,vfdb,tcdb,phi,tmhmm,signalp,secretory
    """

    def __init__(self, work_id):
        super(AnnotationModule, self).__init__(work_id)
        option = [
            {"name": "gene_seq", "type": "infile", "format": "sequence.fasta"},  # 序列
            {"name": "database_list", "type": "string", "default": "cog,kegg,cazy,card,vfdb,tcdb,phi,tmhmm,signalp"},
            # 数据库eggnog,kegg,cazy,card, vfdb, tcdb, phi, tmhmm, signalp
            {"name": "sample", "type": "string"},  # 样品名
            {"name": "kegg", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "gff", "type": "infile", "format": "gene_structure.gff3"}
        ]
        self.add_option(option)
        self.anno_modules = []
        self.anno_dir ={}
        self.cog = self.add_module('bac_comp_genome.anno_cog')
        self.kegg = self.add_module('bac_comp_genome.anno_kegg')
        self.cazy = self.add_module('bac_comp_genome.anno_cazy')
        self.tcdb = self.add_module("bac_comp_genome.anno_tcdb")
        self.card = self.add_module("bac_comp_genome.anno_card")
        self.vfdb = self.add_module("bac_comp_genome.anno_vfdb")
        self.phi = self.add_module("bac_comp_genome.anno_phi")
        self.tmhmm = self.add_tool("bac_comp_genome.tmhmm_anno")
        self.signalp = self.add_module("bacgenome.signalp")
        self.anno_tidy = self.add_tool("bac_comp_genome.anno_tidy")


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("gene_seq").is_set:
            raise OptionError("必须提供基因序列文件")
        if self.option("database_list") == "":
            raise OptionError("必须提供比对的数据名称")
        if self.option("sample") == "":
            raise OptionError("必须提供样品名称")
        return True

    def run_anno(self, run_tool):
        one_run = eval('self.' + run_tool)
        if run_tool == 'signalp':
            query = str(self.option("gene_seq").prop['path'])
        else:
            query = self.option("gene_seq")
        options = {
            "query": query,
            "sample": self.option("sample")
        }
        one_run.set_options(options)
        self.anno_modules.append(one_run)
        self.anno_dir[run_tool] = one_run.output_dir

    def get_all_run(self):
        self.analysis_type = list(set(self.option("database_list").split(',')))
        for i in self.analysis_type:
            self.run_anno(i)
        if len(self.anno_modules) > 1:
            self.on_rely(self.anno_modules, self.run_all_anno)
            self.logger.info(self.anno_modules)
        else:
            self.self.anno_modules[0].on('end', self.set_output)
        for module in self.anno_modules:
            self.logger.info(module)
            module.run()

    def run_all_anno(self):
        opts ={
            "cazy": self.cazy.option("cazy"),
            "gff": self.option("gff"),
            "cog": self.cog.option("cog"),
            "kegg": self.kegg.option("kegg"),
            "card": self.card.option("card"),
            "phi": self.phi.option("phi"),
            "tcdb": self.tcdb.option("tcdb"),
            "vfdb": self.vfdb.option("vfdb"),
            "tmhmm": self.tmhmm.option("tmhmm"),
            "signalp": self.signalp.option("signalp"),
            "signalp1": self.signalp.option("signalp1"),
            "secretory": self.kegg.option("secretory"),
            "sample": self.option("sample"),
        }
        self.anno_tidy.set_options(opts)
        self.anno_tidy.on("end",self.set_output)
        self.anno_tidy.run()

    def set_output(self):
        self.option("kegg", self.kegg.option("kegg"))
        link_dir(self.anno_tidy.output_dir ,self.output_dir + "/" + self.option("sample"))
        link_file(self.vfdb.output_dir + "/" + self.option("sample") + ".vfdb_predict_anno.xls", self.output_dir + "/" + self.option("sample") + "/" + self.option("sample") + ".vfdb_predict_anno.xls")
        link_file(self.vfdb.output_dir + "/" + self.option("sample") + ".vfdb_core_anno.xls",
                  self.output_dir + "/" + self.option("sample") + "/" + self.option("sample") + ".vfdb_core_anno.xls")
        self.end()

    def run(self):
        super(AnnotationModule, self).run()
        self.get_all_run()

    def end(self):
        super(AnnotationModule, self).end()