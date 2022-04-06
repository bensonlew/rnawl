# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2017.04.12
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.rna.merge_kegg_pathway import MergeKeggPathway


class MergeAnnotAgent(Agent):
    """
    将已知（参考基因组）序列和新序列的注释结果合一起
    """
    def __init__(self, parent):
        super(MergeAnnotAgent, self).__init__(parent)
        options = [
            {"name": "gos_dir", "type": "string", "default": None},  # 文件，以；分割
            {"name": "kegg_table_dir", "type": "string", "default": None},
            {"name": "cog_table_dir", "type": "string", "default": None},
            {"name": "pathway_table_dir", "type": "string", "default": None},
            {"name": "database", "type": "string", "default": "go,cog,kegg"},
            {"name": "go2level_out", "type": "outfile", "format": "annotation.go.level2"},
            {"name": "golist_out", "type": "outfile", "format": "annotation.go.go_list"},
            {"name": "kegg_table", "type": "outfile", "format": "annotation.kegg.kegg_table"},
            {"name": "cog_table", "type": "outfile", "format": "annotation.cog.cog_table"}
        ]
        self.add_option(options)
        self.step.add_steps("merge_annot")
        self.on("start", self.step_start)
        self.on("end", self.step_end)
        self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.merge_annot.start()
        self.step.update()

    def step_end(self):
        self.step.merge_annot.finish()
        self.step.update()

    def check_options(self):
        self.database = set(self.option("database").split(","))
        if len(self.database) < 1:
            raise OptionError("至少选择一种注释库")
        for db in self.database:
            if db not in ["go", "cog", "kegg"]:
                raise OptionError("需要合并的注释文件不在支持范围内")
            if db == "go" and not self.option("gos_dir"):
                raise OptionError("缺少go注释的输入文件目录")
            if db == "cog" and not self.option("cog_table_dir"):
                raise OptionError("缺少cog注释table的输入文件目录")
            if db == "kegg":
                if not self.option("kegg_table_dir") and not self.option("pathway_table_dir"):
                    raise OptionError("缺少kegg注释table和pathway_table的输入文件目录")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "注释合并结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["./go2level.xls", "xls", "go注释level2合并文件"],
            ["./gos.list", "xls", "go注释gos合并文件"],
            ["./cog_table.xls", "xls", "cog注释table合并文件"],
            ["./kegg_table.xls", "xls", "kegg注释table合并文件"],
            ["./pathway_table.xls", "xls", "kegg注释pathway合并文件"]
        ])
        super(MergeAnnotAgent, self).end()


class MergeAnnotTool(Tool):
    def __init__(self, config):
        super(MergeAnnotTool, self).__init__(config)
        self._version = '1.0'
        self.database = self.option("database").split(",")
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        self.python = "/miniconda2/bin/python"
        self.merge_scripts = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/merge.py"
        self.goAnnot = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/goAnnot.py"
        self.goSplit = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/goSplit.py"
        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map4.r"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.merge_kegg_pathway = self.config.PACKAGE_DIR + "/rna/merge_kegg_pathway.py"

    def run_merge(self):
        for db in self.database:
            if db == "go":
                tmp = self.option("gos_dir").split(";")
                gos = " ".join(tmp)
                self.merge(dirs=gos, merge_file="query_gos.list", type="go")
                self.option("golist_out", self.work_dir + "/query_gos.list")
                self.run_go_anno()
                self.option("go2level_out", self.work_dir + "/go2level.xls")
                self.logger.info("合并go注释文件完成")
            if db == "cog":
                tmp = self.option("cog_table_dir").split(";")
                cog = " ".join(tmp)
                self.merge(dirs=cog, merge_file="cog_table.xls", type="cog")
                self.option("cog_table", self.work_dir + "/cog_table.xls")
                self.logger.info("合并cog注释文件完成")
            if db == "kegg":
                tmp = self.option("kegg_table_dir").split(";")
                kegg = " ".join(tmp)
                self.merge(dirs=kegg, merge_file="kegg_table.xls", type="kegg")
                self.option("kegg_table", self.work_dir + "/kegg_table.xls")
                self.logger.info("合并kegg注释文件完成")
                r_level_path = self.option("pathway_table_dir").split(";")[0]
                n_level_path = self.option("pathway_table_dir").split(";")[1]
                cmd = "{} {} {} {} {} {} {} {} {}".format(self.python, self.merge_kegg_pathway, self.map_path, self.r_path, self.image_magick, r_level_path, n_level_path, "pathway_table.xls", self.output_dir + "/all_pathways")
                self.logger.info("开始画图")
                self.logger.info(cmd)
                cmd1_obj = self.add_command("merge_png", cmd).run()
                self.wait(cmd1_obj)
                if cmd1_obj.return_code == 0:
                    self.logger.info("画图完成")
                else:
                    self.set_error("画图出错")
                    raise Exception("画图出错")
                # MergeKeggPathway().merge_kegg_pathway(map_path=self.map_path, r_path=self.r_path, image_magick=self.image_magick, r_level_path=r_level_path, n_level_path=n_level_path, all_level_path="pathway_table.xls", all_pathways=self.output_dir + "/all_pathways")
        files = ["go2level.xls", "query_gos.list", "cog_table.xls", "kegg_table.xls", "pathway_table.xls"]
        for f in files:
            if os.path.exists(f):
                linkfile = os.path.join(self.output_dir, f)
                if os.path.exists(linkfile):
                    os.remove(linkfile)
                os.link(f, linkfile)

    def run_go_anno(self):
        cmd1 = "{} {} {} {} {} {}".format(self.python, self.goAnnot, self.option("golist_out").prop["path"], "localhost", self.b2g_user, self.b2g_password)
        self.logger.info("运行goAnnot.py")
        self.logger.info(cmd1)
        cmd1_obj = self.add_command("cmd1", cmd1).run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("运行goAnnot.py完成")
        else:
            cmd1_obj.rerun()
            self.wait(cmd1_obj)
            if cmd1_obj.return_code == 0:
                self.logger.info("运行goAnnot.py完成")
            else:
                self.set_error("运行goAnnot.py出错")
                raise Exception("运行goAnnot.py出错")
        cmd2 = "{} {} {}".format(self.python, self.goSplit, self.work_dir + '/go_detail.xls')
        self.logger.info("运行goSplit.py")
        self.logger.info(cmd2)
        cmd2_obj = self.add_command("cmd2", cmd2).run()
        self.wait(cmd2_obj)
        if cmd2_obj.return_code == 0:
            self.logger.info("运行goSplit.py完成")
        else:
            cmd2_obj.rerun()
            self.wait(cmd2_obj)
            if cmd2_obj.return_code == 0:
                self.logger.info("运行goSplit.py完成")
            else:
                self.set_error("运行goSplit.py出错")
                raise Exception("运行goSplit.py出错")

    def merge(self, dirs, merge_file, type):
        cmd = "{} {} {} {}".format(self.python, self.merge_scripts, dirs, merge_file)
        cmd3_obj = self.add_command("cmd_merge_{}".format(type), cmd).run()
        self.wait(cmd3_obj)
        if cmd3_obj.return_code == 0:
            self.logger.info("文件合并完成")
        else:
            cmd3_obj.rerun()
            self.wait(cmd3_obj)
            if cmd3_obj.return_code == 0:
                self.logger.info("文件合并完成")
            else:
                self.set_error("文件合并未完成")
                raise Exception("文件合并失败")

    def run(self):
        super(MergeAnnotTool, self).run()
        self.run_merge()
        self.end()
