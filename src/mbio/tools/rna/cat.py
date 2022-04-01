# -*- coding: utf-8 -*-
# __author__ = 'shijin'
# modified 2017.04.12
from biocluster.agent import Agent
from biocluster.tool import Tool
import os

class CatAgent(Agent):
    """
    将已知（参考基因组）序列和新序列的注释结果合一起
    """
    def __init__(self, parent):
        super(CatAgent, self).__init__(parent)
        options = [
            {"name": "file1", "type": "string", "default": None},
            {"name": "file2", "type": "string", "default": None},
            {"name": "file3", "type": "outfile", "format": "gene_structure.gtf"}
        ]
        self.add_option(options)
        self.step.add_steps("cat")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.cat.start()
        self.step.update()

    def step_end(self):
        self.step.cat.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 4
        self._memory = "10G"

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "注释合并结果输出目录"],
        # ])
        # result_dir.add_regexp_rules([
        #     ["./go2level.xls", "xls", "go注释level2合并文件"],
        #     ["./gos.list", "xls", "go注释gos合并文件"],
        #     ["./cog_table.xls", "xls", "cog注释table合并文件"],
        #     ["./kegg_table.xls", "xls", "kegg注释table合并文件"]
        # ])
        super(CatAgent, self).end()


class CatTool(Tool):
    def __init__(self, config):
        super(CatTool, self).__init__(config)

    def cat(self):
        cmd = "cat {} {} > ref_new.gtf".format(self.option("file1"), self.option("file2"))
        os.system(cmd)
        self.option("file3", self.work_dir + "/ref_new.gtf")
    
    def run(self):
        super(CatTool, self).run()
        self.cat()
        self.end()
