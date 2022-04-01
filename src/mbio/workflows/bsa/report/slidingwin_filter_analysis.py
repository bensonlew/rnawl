# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.24

import os
import re
import datetime
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SlidingwinFilterAnalysisWorkflow(Workflow):
    """
    交互分析：关联区域定位
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SlidingwinFilterAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "slidingwin_file", "type": "infile", "format": "bsa.vcf"},  # sliding-win.result文件
            {"name": "slid_file", "type": "infile", "format": "bsa.vcf"},  # sliding-win.slid.result文件,用于计算阈值对应关系
            {"name": "region_type", "type": "string", "default": "quantile"},  # 阈值类型，分位值quantile/index值index
            {"name": "region_value", "type": "string"},  # 阈值
            {"name": "i_c_result", "type": "infile", "format": "bsa.vcf"},  # index_calc_result.index
            # {"name": "i_c_result", "type": "string"},  # index_calc_result.index
            {"name": "pop_summary", "type": "infile", "format": "bsa.vcf"},  # pop.summary
            {"name": "p_f_vcf", "type": "infile", "format": "bsa.vcf"},  # pop.final.vcf
            {"name": "wp", "type": "string", "default": ""},  # 野生型亲本名称
            {"name": "mp", "type": "string", "default": ""},  # 突变型亲本名称
            {"name": "wb", "type": "string", "default": ""},  # 野生型混池名称
            {"name": "mb", "type": "string", "default": ""},  # 突变型混池名称
            # {"name": "step", "type": "int"},  # 滑窗策略较小的数值
            {"name": "update_info", "type": "string"},
            {"name": "region_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.slidingwin_filter = self.add_tool("bsa.slidingwin_filter")
        self.region_analysis = self.add_module("bsa.region_analysis")
        self.anno_times = 0

    def check_options(self):
        if not self.option('slidingwin_file'):
            raise OptionError('必须输入sliding-win.result', code="11500301")
        if not self.option('slid_file'):
            raise OptionError('必须输入sliding-win.slid.result', code="11500302")
        if self.option('region_type') not in ["quantile", "index"]:
            raise OptionError('阈值类型只有quantile和index两种', code="11500303")
        if self.option('region_value'):
            if float(self.option('region_value')) < 0 or float(self.option('region_value')) > 1:
                raise OptionError('阈值范围在0-1之间', code="11500304")
        if not self.option('i_c_result'):
            raise OptionError('必须提供index_calc_result结果表', code="11500305")
        if not self.option('pop_summary'):
            raise OptionError('必须提供pop_summary结果表', code="11500306")
        if not self.option('p_f_vcf'):
            raise OptionError('必须提供index_calc_result结果表', code="11500307")
        if not self.option("mb"):
            raise OptionError('必须提供突变型混池mb', code="11500308")
        # if not self.option("step"):
        #     raise OptionError('必须提供滑窗策略较小值')
        return True

    def run_slidingwin_filter(self):
        options = {
            "slidingwin_file": self.option("slidingwin_file"),
            "slid_file": self.option("slid_file"),
            "region_type": self.option("region_type"),
            "region_value": self.option("region_value")
        }
        self.slidingwin_filter.set_options(options)
        self.slidingwin_filter.on("end", self.set_output, "slidingwin_filter")
        self.slidingwin_filter.on("end", self.run_region_analysis)
        self.slidingwin_filter.run()

    def run_region_analysis(self):
        options = {
            "i_c_result": self.option("i_c_result").prop['path'],
            "pop_summary": self.option("pop_summary").prop['path'],
            "p_f_vcf": self.option("p_f_vcf").prop['path'],
            "s_w_select": self.slidingwin_filter.output_dir + "/sliding-win.threshold.select",
            "wp": self.option("wp"),
            "mp": self.option("mp"),
            "wb": self.option("wb"),
            "mb": self.option("mb")
            # "step": self.option("step")
        }
        self.region_analysis.set_options(options)
        self.region_analysis.on("end", self.run_go_anno)
        self.region_analysis.on("end", self.run_kegg_anno)
        self.region_analysis.on("end", self.run_eggnog_anno)
        self.region_analysis.on("end", self.set_output, "region_analysis")
        self.region_analysis.run()

    def run_go_anno(self):
        options = {
            "anno_stat": self.region_analysis.output_dir + "/region_gene/region.threshold.gene.go.stat",
            "anno_type": "go"
        }
        self.go_anno = self.add_tool("bsa.anno_analysis")
        self.go_anno.set_options(options)
        self.go_anno.on("end", self.set_output, "go_anno")
        self.go_anno.run()

    def run_kegg_anno(self):
        options = {
            "anno_stat": self.region_analysis.output_dir + "/region_gene/region.threshold.gene.kegg.stat",
            "anno_type": "kegg"
        }
        self.kegg_anno = self.add_tool("bsa.anno_analysis")
        self.kegg_anno.set_options(options)
        self.kegg_anno.on("end", self.set_output, "kegg_anno")
        self.kegg_anno.run()

    def run_eggnog_anno(self):
        options = {
            "anno_stat": self.region_analysis.output_dir + "/region_gene/region.threshold.gene.eggnog.stat",
            "anno_type": "eggnog"
        }
        self.eggnog_anno = self.add_tool("bsa.anno_analysis")
        self.eggnog_anno.set_options(options)
        self.eggnog_anno.on("end", self.set_output, "eggnog_anno")
        self.eggnog_anno.run()

    def set_output(self, event):
        obj = event["bind_object"]
        if event["data"] == "slidingwin_filter":
            for f in os.listdir(obj.output_dir):
                os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        if event["data"] == "region_analysis":
            for f in os.listdir(obj.output_dir):
                if os.path.isfile(os.path.join(obj.output_dir, f)):
                    os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
                else:
                    if not os.path.exists(os.path.join(self.output_dir, f)):
                        os.mkdir(os.path.join(self.output_dir, f))
                    for f1 in os.listdir(os.path.join(obj.output_dir, f)):
                        os.link(obj.output_dir + "/" + f + "/" + f1, self.output_dir + "/" + f + "/" + f1)
        if re.search(r".+anno", event["data"]):
            self.anno_times += 1
            anno_path = os.path.join(self.output_dir, "region_anno")
            if not os.path.exists(anno_path):
                os.mkdir(os.path.join(anno_path))
            for f in os.listdir(obj.output_dir):
                if os.path.isfile(os.path.join(obj.output_dir, f)):
                    os.link(os.path.join(obj.output_dir, f), os.path.join(anno_path, f))
                else:
                    if not os.path.exists(os.path.join(anno_path, f)):
                        os.mkdir(os.path.join(anno_path, f))
                    for f1 in os.listdir(os.path.join(obj.output_dir, f)):
                        os.link(obj.output_dir + "/" + f + "/" + f1, anno_path + "/" + f + "/" + f1)
        if self.anno_times == 3:
            self.set_db()

    def set_db(self):
        self.logger.info("将数据写到mongo数据库")
        api_region = self.api.api("bsa.bsa_region")
        region_id = self.option("region_id")
        gene_total_path = self.output_dir + "/region_gene/region.threshold.gene.total"
        variant_total_path = self.output_dir + "/region_variant/region.threshold.variant.total"
        vcf_total_path = self.output_dir + "/region_vcf/region.threshold.vcf.total"
        index_quantile_file = self.output_dir + "/quantile.index"
        api_region.add_sg_region_variant(region_id, gene_total_path, variant_total_path)
        api_region.add_sg_region_vcf(region_id, vcf_total_path)
        api_region.add_sg_region_index(region_id, variant_total_path)
        api_region.update_sg_region_quantile(region_id, index_quantile_file)
        api_region_anno = self.api.api("bsa.bsa_region_anno")
        name = "RegionAnnotation_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        anno_id = api_region_anno.add_sg_region_anno(region_id, name)
        gene_total_path = self.output_dir + "/region_gene/region.threshold.gene.total"
        go_stat_detail = self.output_dir + "/region_anno/region.threshold.gene.go.final.stat.detail"
        kegg_stat_detail = self.output_dir + "/region_anno/region.threshold.gene.kegg.final.stat.detail"
        # target_dir = self._sheet.output.strip().split('://')[1]
        target_dir = self._sheet.output.rstrip('/')
        pathway_dir = target_dir + "/region_anno/pathway_dir/"
        eggnog_stat_detail = self.output_dir + "/region_anno/region.threshold.gene.eggnog.final.stat.detail"
        api_region_anno.sg_region_anno_detail(anno_id, gene_total_path)
        api_region_anno.sg_region_anno_go_stat(anno_id, go_stat_detail)
        api_region_anno.sg_region_anno_kegg_stat(anno_id, kegg_stat_detail, pathway_dir)
        api_region_anno.sg_region_anno_eggnog_stat(anno_id, eggnog_stat_detail)
        self.logger.info("导表完成")
        self.end()

    def run(self):
        # self.download_from_s3_()
        self.run_slidingwin_filter()
        super(SlidingwinFilterAnalysisWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SlidingwinFilterAnalysisWorkflow, self).end()
