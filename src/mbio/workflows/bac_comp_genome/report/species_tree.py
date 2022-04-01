# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modifies 20191014


import os,shutil
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from biocluster.file import download, exists
from mbio.packages.bac_comp_genome.common_function import get_fasta, get_sample_from_tree


class SpeciesTreeWorkflow(Workflow):
    """
    物种进化树
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SpeciesTreeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "analysis_type", "type": "string", "default": "s16"},  # s16、ortholog、hous_keeping
            {"name": "type", "type": "string", "default": "port"},  # 核酸nul or 蛋白port
            {"name": "homolog", "type": "infile", "format": "sequence.profile_table"}, #同源蛋白聚类cluster表
            {"name": "gene_path", "type": "string"},
            {"name": "sample_list", "type": "string"},  # 样品信息
            {"name": "16s_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 16s的文件夹
            {"name": "core_gene", "type": "infile", "format": "sequence.fasta"},  # core_gene的合并序列
            {"name": "method", "type": "string", "default": "NJ"},  # NJ or ML
            {"name": "out_group", "type": "string"},  # 外类群在进化文件中名称
            {"name": "bootstrap", "type": "int", "default": 1000},  # 外类群在进化文件中名称
            {"name": "merit", "type": "string", "default": "BIC"},  # 模型选择评估标准，BIC、AIC和AICc
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.species_tree = self.add_module("bac_comp_genome.species_tree")
        self.file_path = self._sheet.output

    def download_files(self):
        self.samples = self.option("sample_list").split(",")
        if self.option("analysis_type") in ['s16']:
            if os.path.exists(self.work_dir + "/16s"):
                shutil.rmtree(self.work_dir + "/16s")
            os.mkdir(self.work_dir + "/16s")
            for sample in self.samples:
                ori_file = self.option("gene_path") + "/" + sample + "/" + sample + "_16S.fna"
                target_file = self.work_dir + "/16s/" + sample + "_16S.ffn"
                download(ori_file, target_file)
            if self.option("out_group"):
                os.link(self.config.SOFTWARE_DIR + "/database/MajorbioDB/" + self.option("out_group") + "/16s/" + self.option("out_group") + "_16S_longest.fa", self.work_dir + "/16s/" +  self.option("out_group") + "_16S.ffn")
        elif self.option("analysis_type") in ['ortholog']:
            if os.path.exists(self.work_dir + "/gene"):
                shutil.rmtree(self.work_dir + "/gene")
            os.mkdir(self.work_dir + "/gene")
            for sample in self.samples:
                if os.path.exists(self.work_dir + "/gene/" + sample):
                    shutil.rmtree(self.work_dir + "/gene/" + sample)
                os.mkdir(self.work_dir + "/gene/" + sample)
                for i in ["_CDS.faa", "_CDS.fna"]:
                    download(self.option("gene_path") + "/" + sample + "/" + sample + i, self.work_dir + "/gene/" + sample + "/" + sample + i)
        elif self.option("analysis_type") in ['house_keeping']:
            get_fasta(self.option("core_gene").prop["path"], self.samples, self.work_dir + "/all.core_gene.fasta")
            if self.option("out_group"):
                os.link(self.config.SOFTWARE_DIR + "/database/MajorbioDB/" + self.option("out_group") + "/house_keeping/housekeeping_gene_31.fasta", self.work_dir + "/" +  self.option("out_group") + ".fa")
                os.system("cat {} >> {}".format(self.work_dir + "/" +  self.option("out_group") + ".fa", self.work_dir + "/all.core_gene.fasta"))

    def run_tree(self):
        if self.option("analysis_type") in ['s16']:
            opts = {
                "analysis_type": self.option("analysis_type"),
                "16s_dir": self.work_dir + "/16s",
                "method": self.option("method"),
                "bootstrap": self.option("bootstrap"),
            }
            if self.option('out_group'):
                opts['out_group'] = self.option('out_group')
        elif self.option("analysis_type") in ['ortholog']:
            opts = {
                "analysis_type": self.option("analysis_type"),
                "type": self.option("type"),
                "homolog": self.option("homolog"),
                "gene_path": self.work_dir + "/gene",
                "sample_list": str(self.option("sample_list")),
                "bootstrap": self.option("bootstrap"),
                "merit": self.option("merit"),
            }
        elif self.option("analysis_type") in ['house_keeping']:
            opts = {
                "analysis_type": self.option("analysis_type"),
                "core_gene": self.work_dir + "/all.core_gene.fasta",
                "bootstrap": self.option("bootstrap"),
                "merit": self.option("merit"),
            }
            if self.option('out_group'):
                opts['out_group'] = self.option('out_group')
        self.species_tree.set_options(opts)
        self.species_tree.on("end", self.set_output)
        self.species_tree.run()


    def run(self):
        self.download_files()
        self.run_tree()
        super(SpeciesTreeWorkflow, self).run()

    def set_output(self):
        self.set_db()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        type = ""
        if self.option("analysis_type") in ['s16']:
            type = "16s"
        elif self.option("analysis_type") in ['ortholog']:
            type = "ortholog"
        elif self.option("analysis_type") in ['house_keeping']:
            type = "house_keeping"
        api_path = self.api.api("bac_comp_genome.common_api")
        tree_id = self.option("main_id")
        if os.path.exists(self.output_dir + "/all." + type + "_" + self.option("method") + ".nwk"):
            os.remove(self.output_dir + "/all." + type + "_" + self.option("method") + ".nwk")
        tree = self.output_dir + '/all.' + type + "_" + self.option("method") + ".nwk"
        file = self.species_tree.output_dir + '/all.tree.nwk'
        os.link(file, tree)
        if os.path.exists(self.output_dir + "/all.assess_result.xls"):
            os.remove(self.output_dir + "/all.assess_result.xls")
        if os.path.exists(self.species_tree.output_dir + "/all.assess_result.xls"):
            os.link(self.species_tree.output_dir + "/all.assess_result.xls", self.output_dir + "/all.assess_result.xls")
        samples = get_sample_from_tree(file)
        path = self.file_path + 'all.' + type + "_" + self.option("method") + ".nwk"
        api_path.add_main_tree(tree, tree_id, update_dic={"type": type, "samples_order": samples, "path": path,"main_id": ObjectId(tree_id), "analysis_type": "species"})
        if os.path.exists(self.output_dir + "/all.assess_result.xls"):
            api_path.add_main_detail(self.output_dir + "/all.assess_result.xls", "module_detail", tree_id, "model,logl,aic,w_aic,aicc,w_aicc,bic,w_bic",main_name="tree_id",convert_dic={"aic": float, "aicc": float, "bic": float})
        self.end()

    def end(self):
        repaths = [
            [".", "", "",0],
        ]
        regexps = [
            [r'.*\.xls$', 'xls', '',0],
            [r'.*\.stat$', 'stat', '',0]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(SpeciesTreeWorkflow, self).end()