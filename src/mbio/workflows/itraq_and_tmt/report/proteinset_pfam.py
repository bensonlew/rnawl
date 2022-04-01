# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import os
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json
import glob


class ProteinsetPfamWorkflow(Workflow):
    """
    蛋白集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetPfamWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "proteinset_files", "type": "string"},
            {"name": "pfam_info", "type": "infile", "format": "itraq_and_tmt.common"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.pfam_stat = self.add_tool("itraq_and_tmt.export_pfam_stat")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/03_Anno/04_Pfam')
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
        super(ProteinsetPfamWorkflow, self).send_log(data)

    def run(self):
        # super(ProteinsetPfamWorkflow, self).run()
        self.pfam_stat.on("end", self.set_db)
        self.run_pfam_stat()
        super(ProteinsetPfamWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_proteinset = self.api.api('itraq_and_tmt.proteinset')
        pfam_file = os.path.join(self.pfam_stat.output_dir, 'pfam_stat.xls')
        try:
            os.link(pfam_file, os.path.join(self.output_dir, 'pfam_stat.xls'))
        except:
            pass

        self.logger.info("开始进行pfam_stat的导表")
        api_proteinset.add_pfam_stat(self.option('main_table_id'), pfam_file)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        pfam_up_down_file = os.path.join(self.output_dir, "pfam_stat.xls")
        chart.chart_pfam_web(pfam_up_down_file)
        chart.to_pdf()

        # move pdf to result dir
        if os.path.exists(os.path.join(self.work_dir, 'pfam.bar.pdf')):
            os.link(os.path.join(self.work_dir, 'pfam.bar.pdf'), self.output_dir + "/pfam.pdf")

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集分析",0],
            ["5_Proteinset/03_Anno/", "", "功能注释", 0],
            ["5_Proteinset/03_Anno/04_Pfam", "", "Pfam注释", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白集Pfam domain分类统计结果目录"],
            ["pfam_stat.xls", " ", "Pfam domain分类统计表"],
            ["pfam.pdf", " ", "Pfam注释统计图"],
        ])
        super(ProteinsetPfamWorkflow, self).end()

    def run_pfam_stat(self):
        opts = {
            "proteinsets": self.option("proteinset_files"),
            "pfam_info": self.option("pfam_info"),
        }
        self.pfam_stat.set_options(opts)
        self.pfam_stat.run()
