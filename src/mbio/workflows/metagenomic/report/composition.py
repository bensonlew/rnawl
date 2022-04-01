# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

"""群落组成分析workflow"""
import os
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import unittest
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
from comm_table import CommTableWorkflow
import json
import shutil

class CompositionWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CompositionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "anno_id", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "params", "type": "string"},
            #{"name": "level_type", "type": "string", "default": ""},
            {"name": "abund_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "others", "type": "float", "default": 0.01},
            {"name": "species_number", "type": "string", "default": "50"},
            {"name": "species_method", "type": "string", "default": ""},
            {"name": "sample_method", "type": "string", "default": ""},
            {"name": "group_method", "type": "string", "default": ""},
            {"name": "normalization", "type": "string", "default": ""},
            {"name": "main_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1},
            {"name": "graphic_type", "type": "string", "default": ""},  # bar,heatmap,circos，bubble
            #{"name": "project_name", "type": "string", "default": ""},  #add by qingchen.zhang@20181204 用于区分与meta的项目
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.composition_analysis = self.add_module("meta.composition.composition_analysis")

    def run_composition_analysis(self):
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            abund_table = self.abundance.option("out_table").prop['path']
        self.logger.info(abund_table)
        self.composition_analysis.set_options({
            "analysis": self.option('graphic_type'),
            "abundtable": abund_table,
            "group": self.option('group'),
            "add_Algorithm":  self.option('group_method'),
            "method":  self.option('species_method'),
            "sample_method":  self.option('sample_method'),
            "species_number":  self.option('species_number'),
            "others":  self.option('others'),
            "level_type": self.option('level_id'),
            "anno_type": self.option('anno_type'),
            "level_color": self.option('level_color'),
            "normalization": self.option('normalization')
        })
        self.composition_analysis.on("end", self.set_db)
        self.composition_analysis.run()

    def replace_0_new_otu_file_path(self):
        """
        将丰度表中的所有丰度为0的在此基础上加上整张丰度表最小值的十分之一,不改变原有结果文件
        :return:
        """
        self.logger.info("正在将结果文件中的丰度值0替换")
        second_list = []
        file_path = self.composition_analysis.output_dir.rstrip('/') + '/heatmap/taxa.table.xls'
        real_zero_otu_new = self.work_dir + "/real_otu_new.xls"
        with open(file_path, 'r') as f, open(real_zero_otu_new, 'w') as w:
            lines = f.readlines()
            w.write("{}".format(lines[0]))
            for line in lines[1:]:
                min_list = []
                line = line.strip().split('\t')
                for i in range(1, len(line)-1):
                    if float(line[i]) != 0.0:
                        line_min = line[i]
                        min_list.append(line_min)
                if min_list:
                    min_line = min(min_list)
                    second_list.append(min_line)
            min_table = float(min(second_list))
            for line in lines[1:]:
                line = line.strip().split('\t')
                data_list = []
                for i in range(1, len(line)-1):
                    if float(line[i]) == 0.0:
                        line[i] = float(line[i]) + float(min_table/10)
                    else:
                        line[i] = float(line[i])
                    data_list.append(line[i])
                w.write("{}\t{}\t{}\n".format(line[0], "\t".join(str(i) for i in data_list), line[-1]))
        return real_zero_otu_new

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        api_composition = self.api.api("metagenomic.composition")
        #api_composition = self.api.composition
        if self.option("graphic_type") in ["heatmap"]:
            specimen_tree = ""
            species_tree = ""
            if self.option("species_method") != "":
                species_tree = self.composition_analysis.output_dir.rstrip('/') + '/heatmap/species_hcluster.tre'
            if self.option("sample_method") != "":
                specimen_tree = self.composition_analysis.output_dir.rstrip('/') + '/heatmap/specimen_hcluster.tre'
            #file_path = self.composition_analysis.output_dir.rstrip('/') + '/heatmap/taxa.table.xls'
            file_path = self.replace_0_new_otu_file_path()
            api_composition.add_composition_detail(file_path, self.option("main_id"), species_tree=species_tree,
                                                   specimen_tree=specimen_tree, group_method=self.option('group_method'), level_color=self.option('level_color'))
        elif self.option("graphic_type") in ["bubble"]:
            file_path = self.composition_analysis.output_dir.rstrip('/') + '/bubble/taxa.percents.table.xls'
            api_composition.add_composition_detail(file_path, self.option("main_id"), species_tree="",
                                                   specimen_tree="", group_method=self.option('group_method'), level_color=self.option('level_color'))
        elif self.option("graphic_type") in ["bar"]:
            file_path = self.composition_analysis.output_dir.rstrip('/') + '/bar/taxa.percents.table.xls'
            api_composition.add_composition_detail(file_path, self.option("main_id"), species_tree="",
                                                   specimen_tree="",group_method=self.option('group_method'))
        elif self.option("graphic_type") in ["circos"]:
            file_path = self.composition_analysis.output_dir.rstrip('/') + '/circos/taxa.percents.table.xls'
            api_composition.add_composition_detail(file_path, self.option("main_id"), species_tree="",
                                                   specimen_tree="", group_method=self.option('group_method'))
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "composition")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = json.loads(self.option('params'))["submit_location"]
            self.logger.info(submit_loc)
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

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))

        if self.option("graphic_type") in ["heatmap"]:
            if os.path.exists(self.output_dir + "/taxa.table.xls"):
                os.remove(self.output_dir + "/taxa.table.xls")
            os.link(self.composition_analysis.output_dir + "/heatmap/taxa.table.xls", self.output_dir + "/taxa.table.xls")
            if os.path.exists(self.composition_analysis.output_dir + "/heatmap/specimen_hcluster.tre"):
                os.link(self.composition_analysis.output_dir + "/heatmap/specimen_hcluster.tre", self.output_dir + "/specimen_hcluster.tre")
            if os.path.exists(self.composition_analysis.output_dir + "/heatmap/species_hcluster.tre"):
                os.link(self.composition_analysis.output_dir + "/heatmap/species_hcluster.tre", self.output_dir + "/species_hcluster.tre")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Heatmap图结果目录", 0, "120099"],
                ["./taxa.table.xls", "xls", "物种/功能丰度结果表", 0, "120168"],
                ["./specimen_hcluster.tre", "tre", "样本聚类树", 0, "120170"],
                ["./species_hcluster.tre", "tre", "物种/功能聚类树", 0, "120171"],
                ["./Heatmap.pdf", "pdf", "群落Heatmap图"]
            ])
        elif self.option("graphic_type") in ["bubble"]:
            if os.path.exists(self.output_dir + "/taxa.table.xls"):
                os.remove(self.output_dir + "/taxa.table.xls")
            os.link(self.composition_analysis.output_dir + "/bubble/taxa.table.xls", self.output_dir + "/taxa.table.xls")
            if os.path.exists(self.output_dir + "/taxa.percents.table.xls"):
                os.remove(self.output_dir + "/taxa.percents.table.xls")
            os.link(self.composition_analysis.output_dir + "/bubble/taxa.percents.table.xls", self.output_dir + "/taxa.percents.table.xls")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Bubble图结果目录", 0, "120281"],
                ["./taxa.table.xls", "xls", "物种/功能丰度结果表", 0, "120168"],
                ["./taxa.percents.table.xls", "xls", "物种/功能相对丰度结果表", 0, "120169"]
            ])
        elif self.option("graphic_type") in ["bar"]:
            if os.path.exists(self.output_dir + "/taxa.table.xls"):
                os.remove(self.output_dir + "/taxa.table.xls")
            os.link(self.composition_analysis.output_dir + "/bar/taxa.table.xls", self.output_dir + "/taxa.table.xls")
            if os.path.exists(self.output_dir + "/taxa.percents.table.xls"):
                os.remove(self.output_dir + "/taxa.percents.table.xls")
            os.link(self.composition_analysis.output_dir + "/bar/taxa.percents.table.xls", self.output_dir + "/taxa.percents.table.xls")
            pdf_dir = ["/nr_pdf", "/cog_pdf", "/kegg_pdf", "/cazy_pdf", "/ardb_pdf", "/card_pdf", "/vfdb_pdf", "/gene_pdf", "/personal_pdf"]
            for dir in pdf_dir:
                if os.path.exists(self.output_dir + dir):
                    os.link(self.output_dir + dir + "/bar.pdf", self.output_dir + "/bar.pdf")
                    shutil.rmtree(self.output_dir + dir)  # 交互分析不存单样本图
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "柱形图结果目录", 0, "120098"],
                ["taxa.table.xls", "xls", "物种/功能丰度结果表", 0, "120168"],
                ["taxa.percents.table.xls", "xls", "物种/功能相对丰度结果表", 0, "120169"],
                ["bar.pdf", "pdf", "群落柱形图"]
            ])
        elif self.option("graphic_type") in ["circos"]:
            if os.path.exists(self.output_dir + "/taxa.table.xls"):
                os.remove(self.output_dir + "/taxa.table.xls")
            if os.path.exists(self.output_dir + "/taxa.percents.table.xls"):
                os.remove(self.output_dir + "/taxa.percents.table.xls")
            os.link(self.composition_analysis.output_dir + "/circos/taxa.table.xls", self.output_dir + "/taxa.table.xls")
            os.link(self.composition_analysis.output_dir + "/circos/taxa.percents.table.xls", self.output_dir + "/taxa.percents.table.xls")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Circos样本与物种或功能关系结果目录", 0, "120100"],
                ["taxa.table.xls", "xls", "物种/功能丰度结果表", 0, "120168"],
                ["taxa.percents.table.xls", "xls", "物种/功能相对丰度结果表", 0, "120169"],
                ["Circos.pdf", "pdf", "Circos图"]
            ])
        super(CompositionWorkflow, self).end()

    def run(self):
        #self.IMPORT_REPORT_DATA = True
        #self.IMPORT_REPORT_DATA_AFTER_END = False
        if self.option("abund_file").is_set:
            self.run_composition_analysis()
        else:
            self.run_abundance(self.run_composition_analysis)
            self.abundance.run()
        super(CompositionWorkflow, self).run()
