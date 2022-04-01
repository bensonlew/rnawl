# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import glob


class WgcnaRelateWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaRelateWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_eigengenes', type='string'),
            dict(name="main_id", type='string'),
            dict(name="trait_path", type='string'),
            dict(name="trait_type", type="string"),
            dict(name="corr_method", type="string"),
            dict(name="block_Rdata", type="infile", format="itraq_and_tmt.common"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("itraq_and_tmt.wgcna.wgcna_relate")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/6_Wgcna')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(WgcnaRelateWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(WgcnaRelateWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("itraq_and_tmt.wgcna")
        # save workflow output path
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        output_dir = self.workflow_output
        dump_tool.update_db_record('sg_wgcna_relate', self.option('main_id'), output_dir=output_dir)
        # add result info
        seq_annot = self.option('exp_eigengenes').split(";")[2]
        dump_tool.add_relate_detail(
            self.tool.work_dir,
            seq_annot,
            self.option('main_id')
        )
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'

        relation_corr = os.path.join(self.tool.work_dir, "module_trait.correlation.xls")
        relation_corr_pvalue = os.path.join(self.tool.work_dir, "module_trait.correlation_pvalues.xls")

        module_stat = os.path.join(self.tool.work_dir, "module_size.stat.xls")

        gene_trait_corr = os.path.join(self.tool.work_dir, "protein_trait.correlation.xls")
        seq_annot = self.option('exp_eigengenes').split(";")[2]
        
        if os.path.exists(relation_corr) and os.path.exists(relation_corr_pvalue) and os.path.exists(module_stat):
            chart.chart_wgcna_relation_corr(relation_corr,relation_corr_pvalue,module_stat)
        if os.path.exists(gene_trait_corr) and os.path.exists(seq_annot) and os.path.exists(relation_corr) and os.path.exists(relation_corr_pvalue):
            chart.chart_wgcna_relation_ms(gene_trait_corr,seq_annot,relation_corr,relation_corr_pvalue)

        chart.to_pdf()

        # move pdf to result dir
        if os.path.exists(os.path.join(self.work_dir, "wgcna.relation_heat.heat_corr.pdf")):
            os.link(os.path.join(self.work_dir, "wgcna.relation_heat.heat_corr.pdf"), os.path.join(self.tool.output_dir, "module_trait_cor.pdf"))
        for pdf_file in glob.glob(self.work_dir + "/wgcna.*.relation_ms.column_conf.pdf"):#chart_wgcna_relation_ms
            file_name = os.path.basename(pdf_file)
            os.link(pdf_file, os.path.join(self.tool.output_dir, file_name[6:-28]+"_MS_bar.pdf"))
        for pdf_file in glob.glob(self.work_dir + "/wgcnarelation_ms.*.scatter.pdf"):#chart_wgcna_relation_ms
            file_name = os.path.basename(pdf_file)
            os.link(pdf_file, os.path.join(self.tool.output_dir, file_name[17:-12]+"_MSGS_plot.pdf"))
        pdf_file = glob.glob(os.path.join(self.tool.output_dir, "*trait.correlation.pdf"))[0]
        if os.path.exists(pdf_file):
            os.link(pdf_file, os.path.join(self.tool.output_dir, "protein_trait_cor.pdf"))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["6_Wgcna", "", "WGCNA",0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "wgcna_relate", 0,  "211250"],
            ["module_trait_cor.pdf", "", "模块与表型相关性热图", 0],
            ["protein_trait_cor.pdf", "", "蛋白与表型相关性热图", 0],
            ["*MS_bar.pdf", "", "MS分析柱图", 0],
            ["*MSGS_plot.pdf", "", "MM-GS散点图", 0],
        ])
        super(WgcnaRelateWorkflow, self).end()

    def run_tool(self):
        options = dict(
            datExpr=self.option('exp_eigengenes').split(";")[0],
            MEs=self.option('exp_eigengenes').split(";")[1],
            traits=self.option('trait_path'),
            corType=self.option('corr_method'),
            block_Rdata=self.option('block_Rdata').prop['path']
        )
        self.tool.set_options(options)
        self.tool.run()
