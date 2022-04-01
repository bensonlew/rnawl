# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.workflow import Workflow
import os
import pandas as pd


class AnnoHmdbWorkflow(Workflow):
    """
    HMDB注释和代谢集HMDB分析共用workflow
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoHmdbWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table,sequence.profile_table"},
            {"name": "metabtable", "type": "infile", "format": "metabolome.metabset,metabolome.metab_abun"},
            {"name": "type", "type": "string", "defaule": "annohmdb"}, # annohmdb, metabsethmbd
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.anno = self.add_tool("metabolome.annotation.anno_hmdb")

    def run(self):
        self.logger.info("start anno hmdb workfolw")
        self.select_anno()
        super(AnnoHmdbWorkflow, self).run()

    def select_anno(self):
        options = {
            "anno_overview": self.option('anno_overview'),
            "metabset": self.option('metabtable'),
            "type" : self.option('type'),
        }
        self.anno.set_options(options)
        self.anno.on('end', self.set_db)
        self.anno.run()

    def set_db(self):
        """
        导表程序
        """
        anno_api = self.api.api('metabolome.anno_hmdb')
        self.output_dir = self.anno.output_dir
        database = self.option('type')
        table = pd.read_table(self.anno.option('level_out').path, sep="\t",header=0)
        if len(table) > 0:
            anno_api.add_level_detail(self.option('main_table_id'), self.anno.option('level_out').path, database =database)
            anno_api.add_stat_detail(self.option('main_table_id'), self.output_dir, database =database)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": self.option("type"),
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("type") == "annohmdb":
            result_dir.add_relpath_rules([
                [".", "", "HMDB化合物分类结果", 0, "150082"],
                ["HmdbLevel.xls", "xls", "HMDB化合物分类层级表", 0, "150083"],
                ["HmdbSuperclass.xls", "xls", "HMDB Superclas统计表", 0, "150084"],
                ["HmdbClass.xls", "xls", "HMDB CLass统计表", 0, "150085"],
                ["HmdbSubclass.xls", "xls", "HMDB Subclass统计表", 0, "150086"]
            ])
        else:
            result_dir.add_relpath_rules([
                [".", "", "代谢集HMDB化合物分类结果", 0, "150087"],
                ["HmdbLevel.xls", "xls", "代谢集HMDB化合物分类层级表", 0, "150088"],
                ["HmdbSuperclass.xls", "xls", "代谢集HMDB Superclas统计表", 0, "150089"],
                ["HmdbClass.xls", "xls", "代谢集HMDB CLass统计表", 0, "150090"],
                ["HmdbSubclass.xls", "xls", "代谢集HMDB Subclass统计表", 0, "150091"]
            ])
        super(AnnoHmdbWorkflow, self).end()
