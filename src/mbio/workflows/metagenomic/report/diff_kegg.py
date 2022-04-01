# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2018/11/28'

from biocluster.workflow import Workflow
import os
from comfun import ComfunWorkflow
# from anno_kegg import AnnoKeggWorkflow
from mbio.packages.metagenomic.common import link_file, link_dir
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class DiffKeggWorkflow(ComfunWorkflow):
    """
    kegg蛋白差异分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffKeggWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "type", "type": "string", "default": "two.side"},
            {"name": "test", "type": "string"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "clean_stat", "type": "infile", "format": "meta.profile"},
            {"name": "method", "type": "string"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.two_group = self.add_tool("metagenomic.diff.diff_kegg")
        self.cal_ppm = self.add_tool("meta.cal_ppm")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run diff kegg workflow")
        if self.option("method") == "ppm" and self.option("clean_stat").is_set:
            self.run_ppm()
        else:
            self.run_filter(self.run_two_group)
        super(DiffKeggWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        link_dir(self.two_group.output_dir + "/pathway_img", self.output_dir + "/pathway_img")
        link_file(self.two_group.option("test_result").path, self.output_dir + "/test_result.xls")
        # api_path = self.api.api('metabolome.ipath')  # edit api path
        # self.logger.info("开始进行导表")
        # api_path.add_ipath_detail(self.option("main_id"), self.ipath.output_dir,
        #                           self.option('metabset').path)  # edit api func
        kegg_api = self.api.api('metagenomic.diff_kegg')
        kegg_api.add_diff_kegg_diff(self.option("main_table_id"), self.output_dir + "/test_result.xls")
        pathway_id = kegg_api.update_main(self.option("main_table_id"), self.filter.recal_abu_tool.output_dir + "/kegg_pathway_profile.xls")  # 必须通过add_diff_kegg_diff获取到丰度最大值
        kegg_api.add_diff_kegg_detail(self.option("main_table_id"), self.output_dir + "/pathway_img/%s.html.mark" % pathway_id)
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_table_id"), "diff_kegg")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_table_id"), "diff_kegg")
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
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
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "代谢通路组间差异分析结果目录",0,"120282"],
            ["pathway_img", "", "通路图片",0,"120283"],
            ["test_result.xls", "", "两组比较检验结果表",0,"120284"],
            ["diffpathway_bar.pdf", "pdf", "代谢通路组间差异检验柱形图"]
        ])
        super(DiffKeggWorkflow, self).end()

    # def run_filter(self):
    #     self.logger.info("start run anno_select module!")
    #     options = {
    #         "gene_anno": self.option("gene_kegg_anno"),
    #         "database": "kegg",
    #         "xml_file": self.option("xml_file"),
    #         "geneset_table": self.option("geneset_table"),
    #         "lowest_level_profile": self.option("lowest_level_profile"),
    #         "gene_profile": self.option("gene_profile"),
    #         "samples": self.option("samples"),
    #         "identity": self.option("identity"),
    #         "align_length": self.option("align_length"),
    #         "level_select": self.option("level_select"),
    #         "abu_num": self.option("abu_num"),
    #         "abu_proportion": self.option("abu_proportion"),
    #         "abu_filter": self.abu_act ,
    #         "fun_filter": self.fun_select,
    #         "gene_filter": self.profile
    #     }
    #     self.filter.set_options(options)
    #     self.filter.on('end', self.run_two_group)
    #     self.filter.run()

    def run_two_group(self):
        options = {
            "test": self.option("test"),
            "kegg_path": self.filter.recal_abu_tool.output_dir,
            "xml_file": self.option("xml_file"),
            "ci": self.option("ci"),
            "group": self.option("group_file"),
            "correction": self.option("correction"),
            "type": self.option("type")
        }
        self.two_group.set_options(options)
        self.two_group.on("end", self.set_db)
        self.two_group.run()

    def run_ppm(self):
        self.logger.info("start calculate ppm>>>>>>>>>>>>>>")
        options = {
            'clean_stat': self.option('clean_stat').prop["path"],
            'geneset_table': self.option('gene_profile').prop["path"],
        }
        self.cal_ppm.set_options(options)
        self.cal_ppm.on("end", self.run_myfilter)
        self.cal_ppm.run()

    def run_myfilter(self):
        self.option("gene_profile", self.cal_ppm.option("out_table").prop["path"])
        self.run_filter(self.run_two_group)
