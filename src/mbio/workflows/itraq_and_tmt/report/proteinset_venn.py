# -*- coding: utf-8 -*-
# __author__ = '刘彬旭'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import time
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import glob
from biocluster.config import Config
from bson.objectid import ObjectId


class ProteinsetVennWorkflow(Workflow):
    """
    蛋白集venn分析
    此workflow 不做任何处理， 交给页面计算， 只用于更新sg_status表
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "proteinset_id", "type": "string"},
            {"name": "type", "type": "string"} # 指示用origin注释还是latest注释
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/02_Venn')
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
        super(ProteinsetVennWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        time.sleep(20)
        self.end()

    def chart(self):
        db = Config().get_mongo_client(mtype="itraq_tmt")[Config().get_mongo_dbname("itraq_tmt")]
        venn_data = []
        for proteinset_id in str(self.option('proteinset_id')).split(','):
            proteinset_id_result1 = db["sg_proteinset"].find_one({"main_id": ObjectId(proteinset_id)})
            proteinset_id_result2= db["sg_proteinset_detail"].find_one({"proteinset_id": ObjectId(proteinset_id)})
            venn_data.append([proteinset_id_result1['name'], proteinset_id_result2['seq_list']])
        # print venn_data
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        chart.chart_proteinset_venn_web(venn_data)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")[0]
        os.link(pdf_file, self.output_dir + "/venn.pdf")

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集分析",0],
            ["5_Proteinset/02_Venn", "", "蛋白venn图", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "", 0],
            ['venn.pdf', 'pdf', 'venn图', 0],
        ])
        super(ProteinsetVennWorkflow, self).end()
