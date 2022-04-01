# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import glob
import os
from shutil import copyfile
import re
from biocluster.core.function import CJsonEncoder
import json
from mbio.packages.denovo_rna_v2.chart import Chart


class ExpVennWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpVennWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="venn_main_id", type='string'),
            dict(name="threshold", type='string'),
            dict(name="group", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("denovo_rna_v2.exp_venn")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/02 Express/05 Exp_Venn')
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
        super(ExpVennWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(ExpVennWorkflow, self).run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        chart.chart_exp_venn(os.path.join(self.tool.output_dir, 'venn_graph.xls'))
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.tool.output_dir + "/" + os.path.basename(p))

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_exp_venn", main_id=self.option('venn_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        # add result info
        graph_table = os.path.join(self.tool.output_dir, 'venn_graph.xls')
        all_exp.add_exp_venn(graph_table, main_id=self.option('venn_main_id'), )
        self.end()

    def end(self):
        self.chart()
        venn_pdf_file = glob.glob(self.work_dir + "/*venn.pdf")[0]
        upset_pdf_file = glob.glob(self.work_dir + "/*upset.pdf")[0]
        if os.path.exists(self.output_dir + "/sample_venn.pdf"):
            os.remove(self.output_dir + "/sample_venn.pdf")
        os.link(venn_pdf_file, self.output_dir + "/sample_venn.pdf")
        if os.path.exists(self.output_dir + "/sample_upset.pdf"):
            os.remove(self.output_dir + "/sample_upset.pdf")
        os.link(upset_pdf_file, self.output_dir + "/sample_upset.pdf")
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["02 Express", "", "表达量分析结果目录", 0],
            ["02 Express/05 Exp_Venn", "", "样本间Venn分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "表达量venn分析结果目录"],
            ['*venn.pdf', 'txt', '样本间Venn图', 0],
            ['*upset.pdf', 'txt', '样本间Upset图', 0],
        ])
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "表达量venn分析结果目录"],
        # ])
        super(ExpVennWorkflow, self).end()

    def run_tool(self):
        options = dict(
            express_matrix=self.option('exp_matrix'),
            group_table=self.option('group'),
            threshold=self.option('threshold'),
        )
        self.tool.set_options(options)
        self.tool.run()
