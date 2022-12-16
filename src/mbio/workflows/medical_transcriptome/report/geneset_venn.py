# -*- coding: utf-8 -*-
# __author__ = 'konghualei, 20170421'
# from mainapp.controllers.project.meta_controller import MetaController
# from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import shutil
import re, os
import glob
from biocluster.workflow import Workflow
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset

class GenesetVennWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/04 GeneSet/05 Venn')
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
        super(GenesetVennWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        self.fire("start")
        # self.set_db()
        self.end()

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录", 0],
            ["04 GeneSet/05 Venn", "", "目标基因集Venn分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "目标基因集Venn分析结果文件", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*venn.pdf', 'txt', '目标基因集Venn图', 0],
            ['*upset.pdf', 'txt', '目标基因集Upset图', 0],
        ])
        super(GenesetVennWorkflow, self).end()

    def set_db(self):
        all_exp = self.api.api("medical_transcriptome.all_exp")
        all_exp.add_venn_tt(
            self.option('main_id'),
        )

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        chart.chart_geneset_venn(self.option("geneset_id"))
        chart.to_pdf()
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/" + os.path.basename(p))