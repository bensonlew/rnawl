# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re
from mbio.packages.labelfree.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import glob
from collections import OrderedDict


class WgcnaPrepareWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaPrepareWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="main_id", type='string'),
            dict(name="me", type='string'),
            dict(name="cv", type="string"),
            dict(name="group_dict", type="string"),
            dict(name="group_id", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("labelfree.wgcna.wgcna_prepare")
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
        super(WgcnaPrepareWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(WgcnaPrepareWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("labelfree.wgcna")
        # add result info
        dump_tool.add_prepare_detail(self.tool.output_dir, main_id=self.option('main_id'), )
        with open(self.work_dir + '/group_info.txt') as f:
            _ = f.readline()
            group_info = dict()
            for line in f:
                s, g = line.strip().split('\t')
                group_info[s] = g
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        exp_matrix =self.workflow_output + '/exp_matrix_after_filtering.txt'
        dump_tool.update_db_record('sg_wgcna_prepare', self.option('main_id'), group_info=group_info, exp_matrix=exp_matrix)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        group_dict = json.loads(self.option('group_dict'), object_pairs_hook=OrderedDict)
        sample_tree = os.path.join(self.tool.output_dir, "sample.cluster.dendrogram.txt")
        wgcna_prepare_curve = os.path.join(self.tool.output_dir, "scale_free_analysis.xls")
        if os.path.exists(sample_tree):
            chart.chart_wgcna_sample_tree(sample_tree,group_dict)
        if os.path.exists(wgcna_prepare_curve):
            chart.chart_wgcna_prepare_curve(wgcna_prepare_curve)
        chart.to_pdf()

        # move pdf to result dir
        for ori_filename,filename in [\
        ["wgcna.sample_tree.heat_corr.pdf","cluster.pdf"],\
        ["wgcna_prepare_average.wgcnascatter.pdf","mean.pdf"],\
        ["wgcna_prepare_adapt.wgcnascatter.pdf","scale.pdf"],\
        ]:
            if os.path.exists(os.path.join(self.work_dir, ori_filename)):
                os.link(os.path.join(self.work_dir, ori_filename), os.path.join(self.tool.output_dir, filename))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["6_Wgcna", "", "WGCNA",0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "wgcna_preprocessing", 0,  "211249"],
            ["cluster.pdf", "", "样本聚类", 0],
            ["scale.pdf", "", "无尺度容适曲线", 0],
            ["mean.pdf", "", "无尺度平均连通度曲线", 0],
        ])
        super(WgcnaPrepareWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
            me=self.option('me'),
            cv=self.option('cv'),
        )
        self.tool.set_options(options)
        self.tool.run()
