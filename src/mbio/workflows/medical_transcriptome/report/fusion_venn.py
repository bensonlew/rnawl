# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart import Chart
import glob


class FusionVennWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FusionVennWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='fusion_matrix', type='infile',format="ref_rna_v2.common"),
            dict(name="group_dict", type="string"),
            dict(name="fusion_venn_main_id", type='string'),
            dict(name="fusion_id", type='string'),
            dict(name="filter_threshold", type='string'),
            dict(name="group", type="string"),
            dict(name="use_group", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("ref_rna_v3.gene_fusion.fusion_venn")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(FusionVennWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("medical_transcriptome.gene_fusion")
        # add result info
        graph_table = os.path.join(self.tool.output_dir, 'fusion_venn_graph.xls')
        all_exp.add_gene_fusion_venn(graph_table, main_id=self.option('fusion_venn_main_id'), )
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        gene_fusion_venn_file =  os.path.join(self.tool.output_dir, 'fusion_venn_graph.xls')
        chart.chart_genefusion_venn(gene_fusion_venn_file)
        chart.to_pdf()
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + os.path.basename(p))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
             [".", "", "基因融合venn分析结果目录"],
         ])
        super(FusionVennWorkflow, self).end()

    def run_tool(self):
        options = dict(
            fusion_matrix=self.option('fusion_matrix'),
            group_table=self.option('group'),
            use_group=self.option('use_group'),
        )
        if self.option("use_group") == "yes":
            options.update({"filter_threshold":self.option('filter_threshold')})
        self.tool.set_options(options)
        self.tool.run()


