# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.workflow import Workflow
import os


class AnnoKeggpWorkflow(Workflow):
    """
    kegg通路注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoKeggpWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_desc"},
            {"name": "organism", "type": "string", "default": "False"},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.anno = self.add_tool("metabolome.annotation.anno_keggp")
        self.output_dir = self.anno.output_dir

    def run(self):
        options = {
            "metab_table": self.option('metab_table'),
            "organism": self.option("organism"),
            "task_id": "_".join(self._sheet.id.split("_")[0:2])
        }
        self.anno.set_options(options)
        self.anno.on('end', self.set_db)
        self.anno.run()
        super(AnnoKeggpWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        anno_api = self.api.api('metabolome.anno_keggp')
        anno_api.add_level_detail(self.option('main_table_id'), self.anno.option('level_out').path)
        anno_api.add_stat_detail(self.option('main_table_id'), self.anno.option('stat_out').path)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "annokeggp",
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
        result_dir.add_relpath_rules([
            [".", "", "KEGG功能通路结果", 0, "150020"],
            ["level.xls", "xls", "Pathway层级表", 0, "150062"],
            ["stat.xls", "xls", "KEGG通路统计表", 0, "150063"],
            ["pathway_img/", "", "KEGG通路图", 0, "150021"]
        ])
        super(AnnoKeggpWorkflow, self).end()
