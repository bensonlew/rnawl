# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.workflow import Workflow
import os


class AnnoKeggcWorkflow(Workflow):
    """
    化合物注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoKeggcWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_desc"},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "database_name", "type": "string", "default": "kegg_compound"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.anno = self.add_tool("metabolome.annotation.anno_keggc")
        self.output_dir = self.anno.output_dir

    def run(self):
        options = {
            "metab_table": self.option('metab_table'),
            "database_name" : self.option("database_name"),
            "task_id": "_".join(self._sheet.id.split("_")[0:2])
        }
        self.anno.set_options(options)
        self.anno.on('end', self.set_db)
        self.anno.run()
        super(AnnoKeggcWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        anno_api = self.api.api('metabolome.anno_keggc')
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
                "submit_loc": "annokeggc",
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
            [".", "", "KEGG化合物分类注释结果", 0, "150017"],
            ["level.xls", "xls", "化合物分类层级表", 0, "150018"],
            ["stat.xls", "xls", "化合物分类统计表", 0, "150019"]
        ])
        super(AnnoKeggcWorkflow, self).end()
