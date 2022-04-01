# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2017.04.13
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess


class KeggUploadAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(KeggUploadAgent, self).__init__(parent)
        options = [
            {"name": "kos_list_upload", "type": "infile", "format": "annotation.upload.anno_upload"},
            {"name": "taxonomy", "type": "string", "default": None},   # kegg数据库物种分类, Animals/Plants/Fungi/Protists/Archaea/Bacteria
            {"name": "kegg_table", "type": "outfile", "format": "annotation.kegg.kegg_table"},
            {"name": "link_bgcolor", "type": "string", "default": "green"},  # 通路图链接官网颜色，约定参考基因组为黄色（yellow），新序列为绿色(green), 两者共有为tomato（红）
            {"name": "png_bgcolor", "type": "string", "default": "#00CD00"}  # 通路图静态图颜色，#00CD00(绿色)，#FFFF00（黄色）
        ]
        self.add_option(options)
        self.step.add_steps('kegg_update')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.kegg_update.start()
        self.step.update()

    def step_end(self):
        self.step.kegg_update.finish()
        self.step.update()

    def check_options(self):
        if not self.option("kos_list_upload").is_set:
            raise OptionError("必须提供kegg注释结果文件", code = "32001201")
        if self.option("taxonomy") not in ["Animals", "Plants", "Fungi", "Protists", "Archaea", "Bacteria", "None", None]:
            raise OptionError("物种类别必须为Animals/Plants/Fungi/Protists/Archaea/Bacteria/None", code = "32001202")

    def set_resource(self):
        self._cpu = 10
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./kegg_table.xls", "xls", "KEGG annotation table"],
            ["./pathway_table.xls", "xls", "Sorted pathway table"],
            ["./kegg_taxonomy.xls", "xls", "KEGG taxonomy summary"]
        ])
        result_dir.add_regexp_rules([
            [r"pathways/ko\d+", 'pdf', '标红pathway图']
        ])
        super(KeggUploadAgent, self).end()


class KeggUploadTool(Tool):

    def __init__(self, config):
        super(KeggUploadTool, self).__init__(config)
        self._version = "2.0"
        self.python = "program/Python/bin/python"
        self.taxonomy_path = self.config.SOFTWARE_DIR + "/database/KEGG/species/{}.ko.txt".format(self.option("taxonomy"))
        # self.kegg_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/kegg_annotation.py"
        self.kegg_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/kegg_annotation_v2.py"
        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map4.r"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"

    def run(self):
        super(KeggUploadTool, self).run()
        self.kegg_annotation()
        self.end()

    def kegg_annotation(self):
        self.logger.info("开始运行kegg注释脚本")
        if not self.option("taxonomy"):
            # taxonomy = self.taxonomy_path
            taxonomy = None
        elif self.option("taxonomy") == "None":
            taxonomy = None
        else:
            taxonomy = self.taxonomy_path
        self.option("kos_list_upload").get_transcript_anno(outdir=self.work_dir + "/kegg.list")
        # self.option("kos_list_upload").get_gene_anno(outdir=self.work_dir + "/gene_kegg.list")
        kegg_ids = self.work_dir + "/kegg.list"
        kegg_table = self.output_dir + '/kegg_table.xls'
        pidpath = self.output_dir + '/pid.txt'
        pathwaydir = self.output_dir + '/pathways'
        pathway_table = self.output_dir + '/pathway_table.xls'
        layerfile = self.output_dir + '/kegg_layer.xls'
        # taxonomyfile = self.output_dir + '/kegg_taxonomy.xls'
        # cmd = "{} {} {} {} {} {} {} {} {} {} {} {}".format(self.python, self.kegg_path, None, kegg_ids, kegg_table, pidpath, pathwaydir, pathway_table, layerfile, taxonomyfile, taxonomy, self.image_magick)
        cmd = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(self.python, self.kegg_path, self.r_path, self.map_path, None, kegg_ids, kegg_table, pidpath, pathwaydir, pathway_table, layerfile, taxonomy, self.option("link_bgcolor"), self.option("png_bgcolor"), self.image_magick)
        command = self.add_command("kegg_anno", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行kegg注释脚本完成")
        else:
            self.set_error("运行kegg注释脚本出错", code = "32001203")
        self.option('kegg_table', self.output_dir + '/kegg_table.xls')
