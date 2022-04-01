# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd
import json
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class ProteinsetIpathWorkflow(Workflow):
    """
    蛋白集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetIpathWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "proteinset_kegg", "type": "string"},
            {"name": "kegg_table", "type": "infile", "format": "labelfree.kegg_table"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "proteinset_id", "type": "string"},
            {"name": "type", "type": "string" , "default": "origin"} # 指示用origin注释还是latest注释
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ipath = self.add_tool("labelfree.proteinset.ipath")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/5_Proteinset/07_Ipath')
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
        super(ProteinsetIpathWorkflow, self).send_log(data)

    def run(self):
        self.ipath.on("end", self.set_db)
        self.run_ipath()
        super(ProteinsetIpathWorkflow, self).run()


    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_proteinset = self.api.api('dia.proteinset')
        self.logger.info("开始进行kegg_class的导表")
        output_file = self.ipath.output_dir + '/gene_ipath_input.xls'
        api_proteinset.add_ipath_detail(self.option("main_table_id"), output_file, self.option('proteinset_kegg'))
        # api_proteinset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        record_id = self.option("main_table_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = api_proteinset.db["sg_proteinset_ipath"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        graph_dir = self.workflow_output
        conn.update({"_id": record_id}, {"$set": {'result_dir': graph_dir}}, upsert=True)
        self.end()

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        result_dir = self.add_upload_dir(self.ipath.output_dir)
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集分析", 0],
            ["5_Proteinset/07_Ipath", "", "Ipath代谢通路图", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白集Ipath结果目录"],
            ["*pdf", "", "iPath代谢通路图"]
        ])
        super(ProteinsetIpathWorkflow, self).end()

    def run_ipath(self):
        opts = {
            "proteinset_kegg": self.option("proteinset_kegg"),
            "kegg_table": self.option("kegg_table"),
        }
        self.ipath.set_options(opts)
        self.ipath.run()
