# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""Venn表计算"""

import os
import datetime
import json
import shutil
import re
from biocluster.workflow import Workflow
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
import gevent
from mbio.packages.meta.save_params import save_params

class VennWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},
            {"name": "update_info", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "samples", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "otu_id", "type": "string"},
            {"name": "venn_id", "type": "string", 'default': ""},
            {"name": "old_venn_id", "type": "string", 'default': ""},## 增加此字段是为了兼容老项目的pie图
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.venn = self.add_tool("graph.venn_table")
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.samples = re.split(',', self.option("samples"))

    def run_venn(self):
        #gevent.sleep(6000)
        options = {
            "otu_table": self.no_zero_otu,
            "group_table": self.option("group_table")
        }
        self.venn.set_options(options)
        self.venn.on('end', self.set_db)
        self.venn.run()

    def run_sort_samples(self):
        """
        运行tool-sort_samples_mg
        功能是为计算pie图计算导表
        :return:
        """
        self.sort_samples.set_options({
            "in_otu_table": self.no_zero_otu,
            "group_table": self.option("group_table"),
            "method": "sum",
        })
        self.sort_samples.on("end", self.run_venn)
        self.sort_samples.run()

    def set_db(self):
        sour = os.path.join(self.venn.work_dir, "output/venn_table.xls")
        dest = os.path.join(self.work_dir, "output")
        shutil.copy2(sour, dest)
        self.logger.info("正在往数据库里插入sg_otu_venn_detail表")
        api_venn = self.api.venn
        # myParams = json.loads(self.sheet.params)
        # name = "venn_table_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        # venn_id = api_venn.create_venn_table(self.sheet.params, myParams["group_id"], self.option("level"), self.option("otu_id"), name)
        venn_id = self.option("venn_id")
        venn_path = os.path.join(self.venn.output_dir, "venn_table.xls")  # modify by zhujuan 20180110 , 2 lines
        venn_graph_path = os.path.join(self.venn.output_dir, "venn_graph.xls")
        asv_table_path = os.path.join(self.sort_samples.output_dir, "taxa.table.xls")
        if self.option("old_venn_id") != "":
            api_venn.add_venn_pie(venn_path, self.option("old_venn_id"),asv_table_path,self.option("group_table").prop['path'])
        else:
            api_venn.add_venn_detail(venn_path, venn_id, self.option("otu_id"), self.option("level"))
            api_venn.add_venn_graph(venn_graph_path, venn_id)
            api_venn.add_venn_pie(venn_path, venn_id,asv_table_path,self.option("group_table").prop['path'])
        # self.add_return_mongo_id("sg_otu_venn", venn_id)
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("venn_id"), "sg_otu_venn")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("venn_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "otu_venn",
                "interaction": 1,
                "main_table": "sg_otu_venn",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        self.no_zero_otu = os.path.join(self.work_dir, "otu.nozero")
        my_sps = self.samples
        self.option("in_otu_table").sub_otu_sample(my_sps, self.no_zero_otu)
        num_lines = sum(1 for line in open(self.no_zero_otu))
        if num_lines < 11:
            self.set_error("Otu表里的OTU数目小于10个！请更换OTU表或者选择更低级别的分类水平！", code="12704301")
        self.run_sort_samples()
        super(VennWorkflow, self).run()

    def end(self):
        self.logger.info("self._sheet.id{}".format(self.id))
        save_params(self.output_dir,self.id)
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Venn图结果目录", 0, "110072"],
            ["venn_table.xls", "xls", "Venn表格", 0, "110073"],
            ["Venn图.pdf", "pdf", "物种Venn图", 0],
        ])
        super(VennWorkflow, self).end()
