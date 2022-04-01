# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.23

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SlidingwinAnalysisWorkflow(Workflow):
    """
    交互分析：标记筛选和分析
    last modified by hongdong @20180611 修改circos的存储路径
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SlidingwinAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pop_vcf", "type": "infile", 'format': "bsa.vcf"},  # pop.final.vcf.gz文件
            {"name": "ref_chrlist", "type": "infile", "format": "bsa.vcf"},  # ref.chrlist文件
            {"name": "pop_summary", "type": "infile", "format": "bsa.vcf"},  # pop.summary文件
            {"name": "wp", "type": "string", "default": ""},  # 野生型亲本名称
            {"name": "mp", "type": "string", "default": ""},  # 突变型亲本名称
            {"name": "wb", "type": "string", "default": ""},  # 野生型混池名称
            {"name": "mb", "type": "string", "default": ""},  # 突变型混池名称
            {"name": "pdep", "type": "int", "default": 10},  # 亲本测序深度，默认10X
            {"name": "bdep", "type": "int", "default": 10},  # 混池测序深度，默认10X
            {"name": "popt", "type": "string", "default": "F2"},  # 群体类型，默认F2（仅qtlseq.pl使用）
            # {"name": "col", "type": "string", "default": "1,2,10"},  # the col of chr pos index
            {"name": "method", "type": "string", "default": "bp"},  # 模式，bp/num
            {"name": "variant_type", "type": "string", "default": "ALL"},  # 只输出变异类型是该类型
            {"name": "win", "type": "string", "default": "2000000"},  # the window size  # modified by hongdong 20180423
            {"name": "step", "type": "string", "default": "10000"},  # the step size
            {"name": "update_info", "type": "string"},
            {"name": "slidingwin_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.slidingwin_analysis = self.add_module("bsa.slidingwin_analysis")
        self.draw_circos = self.add_tool("bsa.draw_circos")
        self.circos_path = ""  # 用于存储circos的画图文件
        self.target_dir = ""

    def check_options(self):
        if not self.option('pop_vcf'):
            raise OptionError('必须输入pop.final.vcf', code="11500201")
        if self.option('method') not in ["bp", "num"]:
            raise OptionError('模式只有bp和num两种', code="11500202")
        if not self.option('mb'):
            raise OptionError('必须输入突变型混池名称', code="11500203")
        return True

    def run_slidingwin_analysis(self):
        options = {
            "pop_vcf": self.option("pop_vcf").prop['path'],
            "wp": self.option("wp"),
            "mp": self.option("mp"),
            "wb": self.option("wb"),
            "mb": self.option("mb"),
            "pdep": self.option("pdep"),
            "bdep": self.option("bdep"),
            "popt": self.option("popt"),
            # "col": self.option("col"),
            "method": self.option("method"),
            "variant_type": self.option("variant_type"),
            "win": self.option("win"),
            "step": self.option("step")
        }
        self.slidingwin_analysis.set_options(options)
        self.slidingwin_analysis.on("end", self.set_output, "slidingwin_analysis")
        self.slidingwin_analysis.on("end", self.run_draw_circos)
        self.slidingwin_analysis.run()

    def run_draw_circos(self):
        options = {
            "p_v_vcf": self.option("pop_vcf").prop['path'],
            "s_w_result": self.slidingwin_analysis.output_dir + "/sliding-win.result",
            "chrlist": self.option("ref_chrlist").prop['path'],
            "pop_summary": self.option("pop_summary").prop['path']
        }
        self.draw_circos.set_options(options)
        self.draw_circos.on("end", self.set_output, "circos")
        self.draw_circos.run()

    def set_output(self, event):
        if event["data"] == "slidingwin_analysis":
            for f in os.listdir(self.slidingwin_analysis.output_dir):
                if os.path.exists(self.output_dir + "/" + f):
                    os.remove(self.output_dir + "/" + f)
                os.link(self.slidingwin_analysis.output_dir + "/" + f, self.output_dir + "/" + f)
        if event["data"] == "circos":
            for f in os.listdir(self.draw_circos.output_dir):
                if os.path.exists(self.output_dir + "/" + f):
                    os.remove(self.output_dir + "/" + f)
                os.link(self.draw_circos.output_dir + "/" + f, self.output_dir + "/" + f)
            self.set_db()

    def set_db(self):
        """这里要修改下后面需要用于计算的文件要上传到远程磁盘"""
        self.logger.info("保存结果到mongo")
        api_slid = self.api.api("bsa.slidingwin")
        stat_file_path = self.output_dir + "/index-calc.result.final.stat"
        slidingwin_id = self.option("slidingwin_id")
        api_slid.add_sg_slidingwin_stat(stat_file_path, slidingwin_id)
        # self.target_dir = self._sheet.output.strip().split('://')[1]
        self.target_dir = self._sheet.output.rstrip("/")
        index_path = self.target_dir + "/index-calc.result.index"
        variant_path = self.target_dir + "/index-calc.result.variant"
        results_path = self.target_dir + "/sliding-win.result"
        slid_path = self.target_dir + "/sliding-win.slid.result"
        api_slid.update_sg_slidingwin_path(slidingwin_id, index_path, variant_path, results_path, slid_path)
        slidingwin = self.api.api("bsa.slidingwin")
        slidingwin.add_sg_circos(slidingwin_id, self.target_dir + "/gene.num.csv",
                                 self.target_dir + "/circos.chrlist", self.output_dir + "/circos.chrlist",
                                 self.target_dir + "/snp.win.csv", self.target_dir + "/indel.win.csv",
                                 self.target_dir + "/sliding.win.csv")
        slidingwin_result_path = self.output_dir + "/sliding-win.result"
        api_slid.add_sg_manhattan(slidingwin_id, slidingwin_result_path, self.option("ref_chrlist"))
        self.logger.info("保存结果到mongo完成")
        self.end()

    def run(self):
        self.run_slidingwin_analysis()
        super(SlidingwinAnalysisWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SlidingwinAnalysisWorkflow, self).end()
