# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '2019.02.11'

from biocluster.workflow import Workflow
import os
import json


class CoreGeneWorkflow(Workflow):
    """
    宏基因组binning模块
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CoreGeneWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "bin_list", "type": "string", "default": ""},  # bin序列序列
            {"name": "genome_list", "type": "string", "default": ""},  # 基因组文件
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.core_bins = []
        self.genome_coregenes = []
        self.all_list = []

    def run(self):
        if self.option("bin_list") != "" and self.option("genome_list") == "":
            self.run_bins()
        elif self.option("bin_list") == "" and self.option("genome_list") != "":
            self.run_genomes()
        elif self.option("bin_list") != "" and self.option("genome_list") != "":
            self.run_bins()
            self.run_genomes()
            self.on_rely(self.all_list, self.set_db)
        super(CoreGeneWorkflow, self).run()

    def run_bins(self):
        samples = eval(self.option("bin_list"))
        self.logger.info(samples)
        for sample,path in samples.items():
            bin_coregene = self.add_module("metagbin.get_bin_coregene")
            opts = ({
                "scaftig": path.replace('\\',""),
                "sample_name": sample
            })
            bin_coregene.set_options(opts)
            self.core_bins.append(bin_coregene)
            if self.option("bin_list") != "" and self.option("genome_list") != "":
                self.all_list.append(bin_coregene)
        if len(self.core_bins) > 1:
            self.on_rely(self.core_bins, self.set_output)
        else:
            self.core_bins[0].on('end', self.set_output)
        for tool in self.core_bins:
            tool.run()

    def run_genomes(self):
        samples = eval(self.option("genome_list"))
        for sample,path in samples.items():
            paths = path.split(',')
            self.logger.info(paths[0].replace('\\',""))
            genome_coregene = self.add_tool("metagbin.get_core_gene")
            opts = ({
                "seq_faa": paths[0].replace('\\',""),
                "seq_gff": paths[1].replace('\\',""),
                "sample_name": sample,
                "method": "genomes"
            })
            genome_coregene.set_options(opts)
            self.genome_coregenes.append(genome_coregene)
            if self.option("bin_list") != "" and self.option("genome_list") != "":
                self.all_list.append(genome_coregene)
        if len(self.genome_coregenes) > 1:
            self.on_rely(self.genome_coregenes, self.set_output)
        else:
            self.genome_coregenes[0].on('end', self.set_output)
        for tool in self.genome_coregenes:
            tool.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        main_id = self.option('main_id')
        self.api_path = self.api.api('metagbin.core_gene')
        self.api_path.add_coregene_detail(main_id, self.output_dir)
        self.end()

    def set_output(self):
        if self.option("bin_list") != "" and self.option("genome_list") == "":
            for tool in self.core_bins:
                files = os.listdir(tool.output_dir)
                for file in files:
                    if os.path.exists(self.output_dir + '/' + file):
                        os.remove(self.output_dir + '/' + file)
                    os.link(tool.output_dir + '/' + file,self.output_dir + '/' + file)
            self.set_db()
        elif self.option("bin_list") == "" and self.option("genome_list") != "":
            for tool in self.genome_coregenes:
                files = os.listdir(tool.output_dir)
                for file in files:
                    if os.path.exists(self.output_dir + '/' + file):
                        os.remove(self.output_dir + '/' + file)
                    os.link(tool.output_dir + '/' + file,self.output_dir + '/' + file)
            self.set_db()
        elif self.option("bin_list") != "" and self.option("genome_list") != "":
            for tool in self.all_list:
                files = os.listdir(tool.output_dir)
                for file in files:
                    if os.path.exists(self.output_dir + '/' + file):
                        os.remove(self.output_dir + '/' + file)
                    os.link(tool.output_dir + '/' + file, self.output_dir + '/' + file)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "CoreGene结果目录"],
            ["*.cor_gene.fa", "", "CoreGene连接合并序列"],
            ["*.result.xls", "", "CoreGene比对看家基因结果表"]
        ])
        super(CoreGeneWorkflow, self).end()