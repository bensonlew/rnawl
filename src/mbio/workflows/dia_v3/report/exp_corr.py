# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
from mbio.packages.dia_v3.chart import Chart
import glob
import os
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re


class ExpCorrWorkflow(Workflow):
    """
    表达量相关性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpCorrWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name='scm', type='string', default='complete'),
            dict(name='scd', type='string', default='correlation'),
            dict(name='corr_method', type='string', default='pearson'),
            dict(name="corr_main_id", type='string'),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("labelfree.exp_corr")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/2_SampleComp/02_SamCorr')
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
        super(ExpCorrWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(ExpCorrWorkflow, self).run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        sample_corr_tree = os.path.join(self.tool.work_dir, 'sample.cluster_tree.txt')
        sample_corr_matrix = os.path.join(self.tool.work_dir, "sample_correlation.xls")
        if os.path.exists(sample_corr_tree) and os.path.exists(sample_corr_matrix):
            chart.chart_exp_corr(sample_corr_matrix, sample_corr_tree)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(os.path.join(self.work_dir, "*.pdf"))
        for p in pdf_file:
            os.link(p, os.path.join(self.tool.output_dir, os.path.basename(p)))

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("dia.all_exp")
        # add result info
        all_exp.add_exp_corr3(self.tool.work_dir, main_id=self.option('corr_main_id'))
        self.end()

    def end(self):
        self.chart()

        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["2_SampleComp", "", "样本比较分析结果目录", 0],
            ["2_SampleComp/02_SamCorr", "", "样本相关性热图结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "相关性分析结果目录"],
            ['./*pdf', 'PDF', '相关性热图'],
        ])
        super(ExpCorrWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
            scm=self.option('scm'),
            scd=self.option('scd'),
            corr_method=self.option('corr_method'),
        )
        self.tool.set_options(options)
        self.tool.run()
