#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == shijin

import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class AnnoUploadModule(Module):
    """
    根据kegg、go的upload_table文件和cog的table文件,进行上传表格文件的注释
    version 1.0
    """
    def __init__(self, work_id):
        super(AnnoUploadModule, self).__init__(work_id)
        options = [
            {"name": "in_cog", "type": "infile", "format": "annotation.cog.cog_table"},
            {"name": "in_kegg", "type": "infile", "format": "annotation.upload.anno_upload"},
            {"name": "in_go", "type": "infile", "format": "annotation.upload.anno_upload"},
            {"name": "database", "type": "string", "default": "go,kegg,string"},
        ]
        self.add_option(options)
        self.kegg = self.add_tool("annotation.kegg_upload")
        self.cog = self.add_tool("annotation.cog_upload")
        self.go = self.add_tool("annotation.go_upload")

    def check_options(self):
        """
        检查参数
        """
        return True

    def run(self):
        super(AnnoUploadModule, self).run()
        self.on_rely([self.go, self.cog, self.kegg], self.end)
        self.run_table()

    def run_table(self):
        if "go" in self.option("database"):
            opts_go = {
                "gos_list_upload": self.option("in_go")
            }
            self.go.set_options(opts_go)
            self.go.run()

        if "kegg" in self.option("database"):
            opts_kegg = {
                "kos_list_upload": self.option("in_kegg")
            }
            self.kegg.set_options(opts_kegg)
            self.kegg.run()
        if "string" in self.option("database"):
            opts_cog = {
                "cogs_list_upload": self.option("in_cog")
            }
            self.cog.set_options(opts_cog)
            self.cog.run()

