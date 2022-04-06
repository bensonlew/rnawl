# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last update liubinxu 20210209
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
import json
from mbio.packages.ref_rna_v2.merge_kegg_pathway import MergeKeggPathway
from mbio.packages.ref_rna_v2.merge import Merge
from mbio.packages.ref_rna_v2.copy_file import CopyFile
import shutil
import pandas as pd


class KeggEditorAgent(Agent):
    """
    kegg 编辑
    """
    def __init__(self, parent):
        super(KeggEditorAgent, self).__init__(parent)
        options = [
            {"name": "input1", "type": "infile", "format": "ref_rna_v2.common"},  # FASTA序列文件
            {"name": "input2", "type": "infile", "format": "ref_rna_v2.common"},  # FASTA序列文件
            {"name": "color_bg1", "type": "string", "default": "#FF0000"},
            {"name": "color_bg2", "type": "string", "default": "#00FF00"},
            {"name": "color_fg11", "type": "string", "default": "#FF0000"},
            {"name": "color_fg12", "type": "string", "default": "#FF0000"},
            {"name": "color_fg21", "type": "string", "default": "#00FF00"},
            {"name": "color_fg22", "type": "string", "default": "#00FF00"},

            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}

        ]
        self.add_option(options)
        self._memory_increase_step = 50
        self.step.add_steps("merge_annot")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.merge_annot.start()
        self.step.update()

    def step_end(self):
        self.step.merge_annot.finish()
        self.step.update()

    def check_options(self):
        pass
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 4
        self._memory = "20G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "注释合并结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["./pathway", "", "kegg编辑图片"]
        ])
        super(KeggEditorAgent, self).end()


class KeggEditorTool(Tool):
    def __init__(self, config):
        super(KeggEditorTool, self).__init__(config)
        self._version = '1.0'
        self.python = "/miniconda2/bin/python"
        self.map_path = self.config.PACKAGE_DIR + "/ref_rna_v2/map4.r"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        if os.path.exists("/usr/bin/convert"):
            self.image_magick = "/usr/bin/convert"
        else:
            self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.kegg_editor = self.config.PACKAGE_DIR + "/tool_lab/kegg_editor.py"
        self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"

    def run_kegg_anno(self):
        # 合并kegg注释
        input1 = self.option("input1").path
        if self.option("input2").is_set:
            input2 = self.option("input2").path
        else:
            input2 = None
        #python ../work/SangerBiocluster2/SangerBiocluster/src/mbio/packages/tool_lab/kegg_editor.py /mnt/ilustre/users/isanger/app/database/KEGG/map_html/ ~/app/program/R-3.3.3/bin/Rscript ~/app/bioinfo/annotation/scripts/map5.r None test1.ko test2.ko \#FFFF00 \#FF00FF \#FF0000,\#0000FF \#00FF00,\#0077FF
        cmd = "{} {} {} {} {} {} {} {} {} {} {} {}".format(
            self.python,
            self.kegg_editor,
            self.html_path,
            self.r_path,
            self.map_path,
            self.image_magick,
            input1,
            input2,
            self.option("color_bg1"),
            self.option("color_bg2"),
            ",".join([self.option("color_fg11"), self.option("color_fg12")]),
            ",".join([self.option("color_fg21"), self.option("color_fg22")]),
        )
        self.logger.info("开始画图")
        self.logger.info(cmd)
        cmd1_obj = self.add_command("editor", cmd, ignore_error=True).run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("画图完成")
        elif cmd1_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("画图出错", code = "33701905")

    def set_output(self):
        for file in os.listdir("pathways"):
            os.link(os.path.join("pathways", file), os.path.join(self.output_dir, file))

    def run(self):
        super(KeggEditorTool, self).run()
        self.run_kegg_anno()
        self.set_output()
        self.end()
