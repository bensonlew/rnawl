# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
from mbio.packages.labelfree.chart import Chart
import glob
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json


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
            dict(name="group", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("labelfree.exp_venn")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/2_SampleComp/03_ExpVenn')
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
        self.run_tool()
        super(ExpVennWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("labelfree.all_exp")
        # add result info
        graph_table = os.path.join(self.tool.output_dir, 'venn_graph.xls')
        all_exp.add_exp_venn(graph_table, main_id=self.option('venn_main_id'), )
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        venn_data = os.path.join(self.tool.output_dir, 'venn_graph.xls')
        if os.path.exists(venn_data):
            chart.chart_exp_venn(venn_data)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(os.path.join(self.work_dir, "*.pdf"))
        for p in pdf_file:
            os.link(p, os.path.join(self.output_dir, os.path.basename(p)))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["2_SampleComp", "", "样本比较分析结果目录", 0],
            ["2_SampleComp/03_ExpVenn", "", "蛋白Venn分析结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "表达量venn分析结果目录"],
            ["./venn_graph.xls", "XLS", "表达量Venn结果详情表", 0],
            ["./*pdf", "PDF", "表达量Venn结果图", 0],
        ])
        super(ExpVennWorkflow, self).end()

    def run_tool(self):
        options = dict(
            express_matrix=self.option('exp_matrix'),
            group_table=self.option('group'),
        )
        self.tool.set_options(options)
        self.tool.run()
