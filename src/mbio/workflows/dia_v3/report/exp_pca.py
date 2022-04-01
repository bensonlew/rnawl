# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import pandas as pd
import os
from mbio.packages.dia_v3.chart import Chart
import glob
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re


class ExpPcaWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpPcaWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="group", type="string"),
            dict(name="pca_main_id", type='string'),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.tool = self.add_tool("itraq_and_tmt.pca")
        self.tool = self.add_tool("labelfree.exp_pca_meta")
        group_dict = json.loads(self.option("group_dict"))
        self.ellipse = None
        i = 0
        for sample_infos in group_dict.values():
            if len(sample_infos) < 3:
                i += 1
            else:
                continue
        if i == 0:
            self.ellipse = self.add_tool("graph.ellipse")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/2_SampleComp/01_SamPca')
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
        super(ExpPcaWorkflow, self).send_log(data)

    def run(self):
        if not self.ellipse is None:
            self.tool.on("end", self.run_ellipse)
            self.ellipse.on("end", self.set_db)
            self.run_tool()
        else:
            self.tool.on("end", self.set_db)
            self.run_tool()
        super(ExpPcaWorkflow, self).run()

    def run_ellipse(self):
        options = {}
        if self.option("group"):
            options['group_table'] = self.option("group")
        options['analysis'] = 'pca'
        options['pc_table'] = self.tool.output_dir + "/pca_sites.xls"
        self.ellipse.set_options(options)
        self.ellipse.run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        sample_pca_sites = os.path.join(self.tool.output_dir, 'pca_sites.xls')
        sample_pca_propor = os.path.join(self.tool.output_dir, 'pca_importance.xls')
        group = os.path.join(self.work_dir, 'group')
        if os.path.exists(sample_pca_sites) and os.path.exists(sample_pca_propor) and os.path.join(group):
            chart.chart_sample_pca(sample_pca_sites, sample_pca_propor, group)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(os.path.join(self.work_dir, "*.pdf"))
        for p in pdf_file:
            if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(p))):
                os.remove(os.path.join(self.tool.output_dir, os.path.basename(p)))
            os.link(p, os.path.join(self.tool.output_dir, os.path.basename(p)))

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("dia.all_exp")
        # add result info
        all_exp.add_exp_pca2(self.tool.output_dir, main_id=self.option('pca_main_id'), )
        if not self.ellipse is None:
            all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(self.option('pca_main_id')))
        else:
            pass
        self.end()

    def end(self):
        # os.remove(self.tool.output_dir + "/pca_rotation.xls")
        # os.rename(self.tool.output_dir + "/pca_rotation_all.xls", self.tool.output_dir + "/pca_rotation.xls")
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["2_SampleComp", "", "样本比较分析结果目录", 0],
            ["2_SampleComp/01_SamPca", "", "PCA分析结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "样本PCA分析结果目录"],
            ["./pca_importance.xls", "xls", "PCA主成分解释表"],
            ["./pca_rotation.xls", "xls", "蛋白相关主成分贡献度"],
            ["./pca_sites.xls", "xls", "PCA分析样本坐标表"],
            ['./*pdf', 'PDF', 'PCA结果图']
        ])
        super(ExpPcaWorkflow, self).end()

    def run_tool(self):
        options = dict(
            # otutable=self.option('exp_matrix'),
            exp_file=self.option('exp_matrix'),
            group_file=self.option('group'),
        )
        self.tool.set_options(options)
        self.tool.run()
