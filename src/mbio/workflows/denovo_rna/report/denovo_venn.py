# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""Venn表计算"""

import os
import shutil
from biocluster.workflow import Workflow


class DenovoVennWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DenovoVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "express_file", "type": "string"},
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},
            {"name": "express_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.venn = self.add_tool("graph.venn_table")
        self.mongodb = self.config.MONGODB + '_rna'

    def run_venn(self):
        options = {
            "otu_table": self.option("express_file").split(',')[1],
            "group_table": self.option("group_table")
        }
        self.logger.info(self.option("express_file"))
        self.venn.set_options(options)
        self.venn.on('end', self.set_db)
        self.venn.run()

    def set_db(self):
        sour = os.path.join(self.venn.work_dir, "output/venn_table.xls")
        dest = os.path.join(self.work_dir, "output")
        shutil.copy2(sour, dest)
        self.logger.info("正在往数据库里插入sg_otu_venn_detail表")
        api_venn = self.api.venn
        api_venn.add_denovo_venn_detail(venn_table=self.output_dir + '/venn_table.xls', venn_id=self.option('main_id'))
        api_venn.add_venn_graph(venn_graph_path=self.venn.work_dir + '/venn_graph.xls', venn_id=self.option('main_id'), project='denovo')
        self.end()

    def run(self):
        self.run_venn()
        super(DenovoVennWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["venn_table.xls", "xls", "Venn表格"]
        ])
        super(DenovoVennWorkflow, self).end()
