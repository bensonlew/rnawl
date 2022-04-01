# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.itraq_and_tmt.chart import Chart
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
            dict(name="level", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("ref_rna_v2.wgcna.wgcna_prepare")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/07 Advanced_Analysis/02 WGCNA')
        self.inter_dirs = []

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(WgcnaPrepareWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="wgcna_prepare", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

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


    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("whole_transcriptome.wgcna")
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
        dump_tool.update_db_record('wgcna_prepare', self.option('main_id'), group_info=group_info, exp_matrix=exp_matrix)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        group_dict = json.loads(self.option('group_dict'), object_pairs_hook=OrderedDict)
        sample_tree = os.path.join(self.tool.output_dir, "sample.cluster.dendrogram.txt")
        wgcna_prepare_curve = os.path.join(self.tool.output_dir, "scale_free_analysis.xls")
        if os.path.exists(sample_tree):
            chart.chart_wgcna_sample_tree(sample_tree, group_dict)
        if os.path.exists(wgcna_prepare_curve):
            chart.chart_wgcna_prepare_curve(wgcna_prepare_curve)
        chart.to_pdf()

        # move pdf to result dir
        for ori_filename, filename in [ \
                ["wgcna.sample_tree.heat_corr.pdf", "cluster.pdf"], \
                ["wgcna_prepare_average.wgcnascatter.pdf", "mean.pdf"], \
                ["wgcna_prepare_adapt.wgcnascatter.pdf", "scale.pdf"], \
                ]:
            if os.path.exists(os.path.join(self.work_dir, ori_filename)):
                os.link(os.path.join(self.work_dir, ori_filename), os.path.join(self.tool.output_dir, filename))
    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录",0],
            ["07 Advanced_Analysis/02 WGCNA", "", "WGCNA分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "WGCNA数据预处理文件", 0],
            ["./powerEstimate_*", "", "soft power 值", 0],
            ["./ignored_gene.list", "", "没有用于聚类的基因列表", 0],
            ["./sample.cluster.*.txt", "", "聚类树形状和分支长度等信息", 0],
            ["./scale_free_analysis.xls", "", "无尺度分析结果", 0],
            ["./exp_matrix_after_filtering.txt", "", "过滤后的表达量表", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
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
