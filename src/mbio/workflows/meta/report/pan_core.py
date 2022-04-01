# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""pan_core OTU计算"""

import shutil
import os
from biocluster.workflow import Workflow
import re
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class PanCoreWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PanCoreWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},
            {"name": "update_info", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "samples", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "main_pan_id", "type": "string"},
            {"name": "main_core_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.pan_core = self.add_tool("meta.otu.pan_core_otu")
        self.samples = re.split(',', self.option("samples"))

    def run_pan_core(self, no_zero_otu):
        if self.option("group_table").prop["is_empty"]:
            options = {
                "in_otu_table": no_zero_otu
            }
        else:
            options = {
                "in_otu_table": no_zero_otu,
                "group_table": self.option("group_table")
            }
        self.pan_core.set_options(options)
        self.pan_core.on('end', self.set_db)
        self.pan_core.run()

    def set_db(self):
        sour = os.path.join(self.pan_core.work_dir, "output/core.richness.xls")
        dest = os.path.join(self.work_dir, "output")
        shutil.copy2(sour, dest)
        sour = os.path.join(self.pan_core.work_dir, "output/pan.richness.xls")
        shutil.copy2(sour, dest)
        self.logger.info("正在写入mongo数据库")
        api_pan_core = self.api.pan_core
        pan_path = self.pan_core.option("pan_otu_table").prop['path']
        core_path = self.pan_core.option("core_otu_table").prop['path']
        api_pan_core.add_pan_core_detail(pan_path, self.option('main_pan_id'))
        api_pan_core.add_pan_core_detail(core_path, self.option('main_core_id'))
        #self.end()
        self.save_pdf1()

    def save_pdf1(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_pan_id"), "sg_otu_pan_core")
            self.figsave1 = self.add_tool("meta.fig_save")
            self.figsave1.on('end', self.save_pdf2)
            self.figsave1.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_pan_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "otu_pan_core",
                "interaction": 1,
                "main_table": "sg_otu_pan_core",
            })
            self.figsave1.run()
        else:
            self.end()

    def save_pdf2(self):
        name = get_name(self.option("main_core_id"), "sg_otu_pan_core")
        self.figsave2 = self.add_tool("meta.fig_save")
        self.figsave2.on('end', self.end)
        self.figsave2.set_options({
            "task_id": "_".join(self.sheet.id.split('_')[:2]),
            "table_id": self.option("main_core_id"),
            "table_name": name,
            "project": "meta",
            "submit_loc": "otu_pan_core",
            "interaction": 1,
            "main_table": "sg_otu_pan_core",
        })
        self.figsave2.run()

    def run(self):
        no_zero_otu = os.path.join(self.work_dir, "otu.nozero")
        my_sps = self.samples
        self.option("in_otu_table").sub_otu_sample(my_sps, no_zero_otu)
        num_lines = sum(1 for line in open(no_zero_otu))
        if num_lines < 11:
            self.set_error("OTU表里的OTU数目小于10个！请更换OTU表或者选择更低级别的分类水平！", code="12703201")   #将Otu改成了OTU modified by hongdongxuan 20170310
        self.run_pan_core(no_zero_otu)
        super(PanCoreWorkflow, self).run()

    def end(self):
        if self.pdf_status:
            if os.path.exists(self.output_dir + "Pan曲线图.pdf"):
                os.remove(self.output_dir + "/Pan曲线图.pdf")
            if os.path.exists(self.figsave1.output_dir + "/1曲线图.pdf"):
                os.link(self.figsave1.output_dir + "/1曲线图.pdf", self.output_dir + "/" + "Pan曲线图.pdf")
            if os.path.exists(self.output_dir + "/Core曲线图.pdf"):
                os.remove(self.output_dir + "/Core曲线图.pdf")
            if os.path.exists(self.figsave2.output_dir + "/2曲线图.pdf"):
                os.link(self.figsave2.output_dir + "/2曲线图.pdf", self.output_dir + "/Core曲线图.pdf")
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Pan/Core结果目录", 0, "110032"],
            ["core.richness.xls", "xls", "Core 表格", 0, "110034"],
            ["pan.richness.xls", "xls", "Pan 表格", 0, "110033"],
            ["Pan曲线图.pdf", "pdf", "Pan曲线图", 0, ""],
            ["Core曲线图.pdf", "pdf", "Core曲线图", 0, ""]
        ])
        print self.get_upload_files()
        super(PanCoreWorkflow, self).end()
