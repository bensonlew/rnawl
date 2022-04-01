# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import os
from mbio.packages.labelfree.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json
import glob


class ProteinsetSublocWorkflow(Workflow):
    """
    蛋白集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetSublocWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "proteinset_files", "type": "string"},
            {"name": "subloc_info", "type": "infile", "format": "labelfree.common"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.subloc_stat = self.add_tool("labelfree.export_subloc_stat")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/03_Anno/05_SubLoc')
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
        super(ProteinsetSublocWorkflow, self).send_log(data)

    def run(self):
        # super(ProteinsetSublocWorkflow, self).run()
        self.subloc_stat.on("end", self.set_db)
        self.run_subloc_stat()
        super(ProteinsetSublocWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_proteinset = self.api.api('labelfree.proteinset')
        subloc_file = os.path.join(self.subloc_stat.output_dir, 'subloc_stat.xls')
        try:
            os.link(subloc_file, os.path.join(self.output_dir, 'subloc_stat.xls'))
        except:
            pass

        self.logger.info("开始进行subloc_stat的导表")
        api_proteinset.add_subloc_stat(self.option('main_table_id'), subloc_file)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        subloc_up_down_file = os.path.join(self.output_dir, "subloc_stat.xls")
        chart.chart_subloc_web(subloc_up_down_file)
        chart.to_pdf()

        # move pdf to result dir
        if os.path.exists(os.path.join(self.work_dir, 'subloc.bar.pdf')):
            os.link(os.path.join(self.work_dir, 'subloc.bar.pdf'), self.output_dir + "/subloc.pdf")

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集分析",0],
            ["5_Proteinset/03_Anno", "", "功能注释", 0],
            ["5_Proteinset/03_Anno/05_SubLoc", "", "Subloc注释", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白集亚细胞定位分类统计结果目录"],
            ["subloc_stat.xls", " ", "亚细胞定位分类统计表"],
            ["subloc.pdf", " ", "亚细胞定位分类统计图"],
        ])
        super(ProteinsetSublocWorkflow, self).end()

    def run_subloc_stat(self):
        opts = {
            "proteinsets": self.option("proteinset_files"),
            "subloc_info": self.option("subloc_info"),
        }
        self.subloc_stat.set_options(opts)
        self.subloc_stat.run()
