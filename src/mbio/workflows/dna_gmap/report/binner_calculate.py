# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180718

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class BinnerCalculateWorkflow(Workflow):
    """
    BinMarker分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BinnerCalculateWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pop_final_vcf", "type": "infile", "format": "bsa.vcf"},  # pop.final.vcf.gz文件，用于染色体bin统计
            {"name": "genotype_matrix", "type": "infile", "format": "dna_gmap.marker"},  # 分型矩阵，此处为pop.filtered.marker,遗传标记筛选的结果文件
            {"name": "marker_upload", "type": "infile", "format": "dna_gmap.marker"},  # 客户上传的标记列表
            {"name": "pop_type", "type": "string"},  # 群体类型，F2，F1，BC，CP，RIL  有CP是F1，无CP是F2,
            {"name": "window_size", "type": "float", "default": 2000},  # window size(kb)
            {"name": "window_step", "type": "float", "default": 100},  # step size(kb)
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("pop_final_vcf").is_set:
            raise OptionError("请设置pop.final.vcf文件", code="14800201")
        if not self.option("genotype_matrix").is_set:
            raise OptionError("请设置分型矩阵genotype_matrix", code="14800202")
        if not self.option("pop_type"):
            raise OptionError("请设置群体类型", code="14800203")
        if self.option("pop_type") not in ["F1", "F2", "BC", "CP", "RIL"]:
            raise OptionError("群体类型只能在F1/F2/BC/CP/RIL内", code="14800204")

    def run_binner_calculate(self):
        options = {
            "pop_final_vcf": self.option("pop_final_vcf"),
            "genotype_matrix": self.option("genotype_matrix"),
            "pop_type": self.option("pop_type"),
            "window_size": self.option("window_size"),
            "window_step": self.option("window_step")
        }
        if self.option("marker_upload"):
            options["marker_upload"] = self.option("marker_upload")
        self.binner_calculate = self.add_tool("dna_gmap.binner_calculate")
        self.binner_calculate.set_options(options)
        self.binner_calculate.on("end", self.set_output, "binmarker")
        self.binner_calculate.run()

    def set_output(self, event):
        obj = event['bind_object']
        for f in os.listdir(obj.output_dir):
            old = os.path.join(obj.output_dir, f)
            new = os.path.join(self.output_dir, f)
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)
        self.set_db()

    def set_db(self):
        """
        将结果导入mongo数据库
        """
        self.logger.info("将结果导入mongo数据库")
        binner_api = self.api.api("dna_gmap.binner_calculate")
        binmarker_id = self.option("main_id")
        bin_stat = os.path.join(self.output_dir, "bin_stat.xls")
        bin_info = os.path.join(self.output_dir, "bin_info.xls")
        bin_pos = os.path.join(self.output_dir, "bin_pos.xls")
        filter_marker_path = os.path.join(self._sheet.output.rstrip('/'), "Total.bin.marker")
        binner_api.add_sg_binmarker_bin(binmarker_id, bin_stat)
        binner_api.add_sg_binmarker_var(binmarker_id, bin_info)
        binner_api.add_sg_binmarker_var_detail(binmarker_id, bin_pos)
        binner_api.update_sg_binmarker_path(binmarker_id, filter_marker_path)
        self.end()

    def run(self):
        self.run_binner_calculate()
        super(BinnerCalculateWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BinnerCalculateWorkflow, self).end()
