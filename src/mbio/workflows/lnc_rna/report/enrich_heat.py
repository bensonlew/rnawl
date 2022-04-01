# -*- coding: utf-8 -*-
# 2019-03-21
import json
import time

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class EnrichHeatWorkflow(Workflow):
    """
    基因集功能分类分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnrichHeatWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "data", "type": "string"},
            {"name": "geneset_type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/04 GeneSet/09 Heatmap')
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
        super(EnrichHeatWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_run_log()
        self.set_db()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_geneset_go_enrich_heatmap", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        self.logger.debug('=============== insert data to table ===============')
        api_geneset = self.api.api('lnc_rna.lnc_rna_geneset')
        self.logger.info("开始进行kegg_class的导表")
        data_file = self.option('data')
        # main_table_id, data_path, table_name
        update_info = json.loads(self.option('update_info'))
        main_table_id = self.option('main_table_id')
        main_table = update_info[main_table_id]
        api_geneset.add_enrich_heatmap(main_table_id, data_file, main_table)
        self.logger.debug('========= insert data to table completed ===========')
        conn = api_geneset.db[main_table]
        conn.update({"_id": ObjectId(main_table_id)}, {"$set": {'status': 'end'}}, upsert=True)
        self.end()

    def end(self):
        self.logger.debug('start end function')
        data_file = self.option('data')
        new_file = os.path.join(self.output_dir, os.path.basename(data_file))
        os.system('ln {} {}'.format(data_file, new_file))
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录", 0],
            ["04 GeneSet/09 Heatmap", "", "富集热图", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "富集热图文件", 0, ],
            ['./heatmap_data.xls', '', '富集热图结果表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        self.logger.debug('end function')
        super(EnrichHeatWorkflow, self).end()
