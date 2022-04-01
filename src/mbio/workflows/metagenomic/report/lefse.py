# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""lefse分析"""

from biocluster.workflow import Workflow
import os
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class LefseWorkflow(CommTableWorkflow):
    """
    报告中调用lefse分析时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LefseWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "lefse_input", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "anno_table", "type": "infile", "format": "meta.profile"},
            {"name": "params", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "anno_id", "type": "string"},
            {"name": "level_id", "type": "string"},
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},
            #{"name": "func_table", "type": "infile", "format": "meta.profile"},
            #{"name": "tax_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            {"name": "level_type", "type": "string"},
            {"name": "level_type_name", "type": "string", "default": ""},
            {"name": "lowest_level", "type": "string", "default": ""},
            {"name": "lefse_gname", "type": "string"},
            {"name": "second_group_detail", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "lda_filter", "type": "float", "default": 0},
            {"name": "strict", "type": "int", "default": 0},
            {"name": "group_id", "type": "string"},
            {"name": "lefse_type", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "start_level", "type": "int", "default": 3},
            {"name": "end_level", "type": "int", "default": 7},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            {"name":"domain_type","type":"string","default": "Bacteria"},
            {"name": "clean_stat", "type": "infile", "format": "meta.profile"},
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
            {"name": "method", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.lefse = self.add_module("meta.diff.lefse")

    def run_lefse(self):
        if self.option("method") == "ppm" and self.option("clean_stat").is_set and not self.option("anno_type") == "go":
            geneset_table = self.cal_ppm.option("out_table").prop["path"]
        elif self.option("method") == "ppm" and self.option("clean_stat").is_set and self.option("anno_type") == "go":
            geneset_table = self.abundance.option("out_table").prop["path"]
        else:
            geneset_table = self.option('geneset_table')
        if not self.option('lefse_input').is_set:
            if self.option('lefse_type') in ['metagenome_taxon']:
                options = {
                 #  "tax_table": self.option("tax_table"),
                   "lefse_group": self.option("group"),
                   "lda_filter": self.option("lda_filter"),
                   "strict": self.option("strict"),
                   "lefse_type":self.option('lefse_type'),
                   "lefse_gname": self.option("lefse_gname"),
                   "start_level": self.option("start_level"),
                   "end_level": self.option("end_level"),
                   "geneset_table": geneset_table,
                   "anno_table": self.option("anno_table"),
                   "domain_type": self.option("domain_type"),
                }
            elif self.option('lefse_type') in ['metagenome_func']:
                if self.option("anno_type") == "go":
                    input = self.abundance.option("out_table").prop["path"]
                    options ={
                        "lefse_input": input,
                        "lefse_group": self.option("group"),
                        "lda_filter": self.option("lda_filter"),
                        "strict": self.option("strict"),
                        "lefse_type": self.option('lefse_type'),
                        "lefse_gname": self.option("lefse_gname"),
                        "start_level": self.option("start_level"),
                        "end_level": self.option("end_level"),
                        'geneset_table': geneset_table
                    }
                elif not self.option('gene_list').is_set:
                    options = {
                       "lefse_group": self.option("group"),
                       "lda_filter": self.option("lda_filter"),
                       "lefse_type": self.option('lefse_type'),
                       "strict": self.option("strict"),
                       "lefse_gname": self.option("lefse_gname"),
                       'anno_table': self.option("anno_table"),
                       'geneset_table': geneset_table,
                       'gene_list': self.option('gene_list'),
                       'level_type': self.option('level_type'),
                       'level_type_name': self.option('level_type_name'),
                       'lowest_level': self.option('lowest_level'),
                }
                else:
                    options = {
                        "lefse_group": self.option("group"),
                        "lda_filter": self.option("lda_filter"),
                        "lefse_type": self.option('lefse_type'),
                        "strict": self.option("strict"),
                        "lefse_gname": self.option("lefse_gname"),
                        'geneset_table': geneset_table,
                        'gene_list': self.option('gene_list'),
                    }
        else:
            options ={
                "lefse_input": self.option("lefse_input"),
                "lefse_group": self.option("group"),
                "lda_filter": self.option("lda_filter"),
                "strict": self.option("strict"),
                "lefse_type": self.option('lefse_type'),
                "lefse_gname": self.option("lefse_gname"),
                "start_level": self.option("start_level"),
                "end_level": self.option("end_level"),
            }
        self.lefse.set_options(options)
        self.lefse.on("end", self.set_db)
        self.lefse.run()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.lefse.output_dir))
        result_dir = self.add_upload_dir(self.lefse.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "LEfSe分析结果目录", 0, "120192"],
            # ["./lefse_LDA.cladogram.png", "png", "LEfSe分析cladogram结果图片"],
            # ["./lefse_LDA.png", "png", "LEfSe分析LDA图片"],
            ["./lefse_LDA.xls", "xls", "LEfSe分析lda数据表", 0, "120193"],
            ["lefse_cladogram.pdf", "pdf", "多级物种层级树图"],
            ["lefse_LDA.pdf", "pdf", "LDA判别柱形图"],
        ])
        super(LefseWorkflow, self).end()

    def set_db(self):
        """
        保存两组比较分析的结果表保存到mongo数据库中
        """
        if self.option('lefse_type') == 'meta_func':
            api_lefse = self.api.api('stat_test')
            lefse_path = self.lefse.output_dir + '/lefse_LDA.xls'
            if not os.path.isfile(lefse_path):
                raise Exception("找不到报告文件:{}".format(lefse_path))
            api_lefse.add_species_difference_lefse_detail(file_path=lefse_path, table_id=self.option("main_id"))
        else:
            api_lefse = self.api.api('metagenomic.lefse')
            lefse_path = self.lefse.output_dir + '/lefse_LDA.xls'
            if not os.path.isfile(lefse_path):
                raise Exception("找不到报告文件:{}".format(lefse_path))
            api_lefse.add_lefse_detail(file_path=lefse_path, lefse_id=self.option("main_id"))
        #if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "lefse")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = json.loads(self.option('params'))["submit_location"]
            # self.logger.info(submit_loc)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        if self.option("method") == "ppm" and self.option("clean_stat").is_set:
            if self.option("anno_type") == "go":
                self.run_abundance(self.run_lefse)
                self.abundance.run()
            else:
                self.run_cal_ppm(rely=self.run_lefse, is_run=True)
        else:
            if self.option("anno_type") == "go":
                self.run_go_abu(rely=self.run_lefse, is_run=True, level=self.option("level_type"))
            else:
                self.run_lefse()
        super(LefseWorkflow, self).run()