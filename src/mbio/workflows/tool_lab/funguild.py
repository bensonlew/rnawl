# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_file,link_dir
import os

class FunguildWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FunguildWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_table", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.funguild = self.add_tool('tool_lab.funguild')

    def check_options(self):
        if not self.option("data_table").is_set:
            raise OptionError("请传入数据表！", code="")

    def run(self):
        self.run_funguild()
        super(FunguildWorkflow, self).run()

    def run_funguild(self):
        tax_abund_table = self.option("data_table").prop['path']
        self.funguild.set_options({
            "taxon_table": tax_abund_table,
        })
        self.funguild.on("end", self.set_db)
        self.funguild.run()

    def set_db(self):
        link_dir(self.funguild.output_dir, self.output_dir)
        api_main = self.api.api("tool_lab.common")
        api_main.add_main_table("funguild", main_id = self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "FUNGuild功能预测结果目录", 0],
            ["./Funguild.txt", "xls", "FUNGuild功能预测结果表", 0],
            ["./FUNGuild_guild.txt", "txt", "Guild功能分类统计表", 0],
            ["./FUNGuild_Trophic_mode.txt", "txt", "Trophic_mode丰度表", 0]
        ])
        super(FunguildWorkflow, self).end()