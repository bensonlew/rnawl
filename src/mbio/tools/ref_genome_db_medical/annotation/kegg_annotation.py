# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan'
# modified 2016.11.28
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import re
from biocluster.core.exceptions import OptionError
import xml.etree.ElementTree as ET
import subprocess


class KeggAnnotationAgent(Agent):
    """
    to perform KEGG annotation
    author:chenyanyan
    modified at 20161128
    """

    def __init__(self, parent):
        super(KeggAnnotationAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "taxonomy", "type": "string", "default": None},   # kegg数据库物种分类, Animals/Plants/Fungi/Protists/Archaea/Bacteria
            {"name": "kegg_table", "type": "outfile", "format": "ref_rna_v2.kegg_table"},
            {"name": "known_ko", "type": "string", "default": None},
            {"name": "kegg_species", "type": "string", "default": None},
            {"name": "link_bgcolor", "type": "string", "default": "green"},  # 通路图链接官网颜色，约定参考基因组为黄色（yellow），新序列为绿色(green), 两者共有为tomato（红）
            {"name": "png_bgcolor", "type": "string", "default": "#00CD00"}  # 通路图静态图颜色，#00CD00(绿色)，#FFFF00（黄色）
        ]
        self.add_option(options)
        self._memory_increase_step = 20
        self.step.add_steps('kegg_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.kegg_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.kegg_annotation.finish()
        self.step.update()

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("必须提供BLAST结果文件", code = "33701501")
        if self.option("taxonomy") not in ["Animals", "Plants", "Fungi", "Protists", "Archaea", "Bacteria", "None", None, "All"]:
            raise OptionError("物种类别必须为Animals/Plants/Fungi/Protists/Archaea/Bacteria/None", code = "33701502")

    def set_resource(self):
        self._cpu = 2
        file_size = float(os.path.getsize(self.option('blastout').prop['path'])) / 1024 / 1024
        mem = int(float(file_size)/1024 * 18) + 2
        mem = min(mem, 180)
        self._memory = '{}G'.format(mem + 16)

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
        super(KeggAnnotationAgent, self).end()


class KeggAnnotationTool(Tool):

    def __init__(self, config):
        super(KeggAnnotationTool, self).__init__(config)
        self._version = "2.0"
        self.python = "miniconda2/bin/python"
        self.taxonomy_path = self.config.SOFTWARE_DIR + "/database/KEGG/species/{}.ko.txt".format(self.option("taxonomy"))
        # self.kegg_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/kegg_annotation.py"
        #self.config.PACKAGE_DIR + "/denovo_rna_v2/new_annotation_query.py"
        self.kegg_path3 = self.config.PACKAGE_DIR + "/ref_rna_v2/kegg_annotation_v3.py"
        self.kegg_path2 = self.config.PACKAGE_DIR + "/ref_genome_db_v2/kegg_annotation_v2.py"
        self.map_path = self.config.PACKAGE_DIR + "/ref_rna_v2/map4.r"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"

        if os.path.exists("/usr/bin/convert"):
            self.image_magick = "/usr/bin/convert"
        else:
            self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.3/lib64/R/lib:$LD_LIBRARY_PATH"
        self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.3/lib64/R/"
        self.set_environ(R_LIBS=self.config.SOFTWARE_DIR + '/program/R-3.3.3/lib64/R/lib')
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run(self):
        super(KeggAnnotationTool, self).run()
        self.kegg_annotation()
        self.end()

    def kegg_annotation(self):
        self.logger.info("运行kegg注释脚本")
        if not self.option("taxonomy"):
            taxonomy = None
        elif self.option("taxonomy") == "None":
            taxonomy = None
        else:
            taxonomy = self.taxonomy_path
        blast_xml = self.option('blastout').prop['path']
        kegg_table = self.output_dir + '/kegg_table.xls'
        pidpath = self.output_dir + '/pid.txt'
        pathwaydir = self.output_dir + '/pathways'
        image_magick = self.image_magick
        pathway_table = self.output_dir + '/pathway_table.xls'
        layerfile = self.output_dir + '/kegg_layer.xls'
        # taxonomyfile = self.output_dir + '/kegg_taxonomy.xls'
        if self.option('kegg_species'):
            genome_path = os.path.join(self.config.SOFTWARE_DIR, "database/KEGG/genome2.xls")
            genome_id = ''
            genome_abr = ''
            genome = ''
            species_path_dir = ''
            with open (genome_path, 'r' ) as f:
                lines = f.readlines()
                for line in lines:
                    genome_id = re.sub("gn:", "", line.split("\t")[0])
                    genome_abr = line.split("\t")[1].split(',')[0]
                    genome = line.split("\t")[1].split(';')[-1].strip()
                    if self.option("kegg_species") == genome_abr or self.option("kegg_species") == genome:
                        break
            if genome_id:
                species_path_dir = os.path.join(self.config.SOFTWARE_DIR, "database/KEGG/kegg_2017-05-01/kegg/pathway/organisms", genome_abr )
            else:
                self.logger.info("未找到相关的物种信息")
            cmd = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(self.python, self.kegg_path3, self.r_path, self.map_path, blast_xml, None, kegg_table, pidpath, pathwaydir, pathway_table, layerfile, taxonomy, self.option("link_bgcolor"), self.option("png_bgcolor"), self.image_magick, genome_abr, species_path_dir, self.html_path, self.option("known_ko"))
        else:
            cmd = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(self.python, self.kegg_path2, self.r_path, self.map_path, blast_xml, None, kegg_table, pidpath, pathwaydir, pathway_table, layerfile, taxonomy, self.option("link_bgcolor"), self.option("png_bgcolor"), self.image_magick, self.html_path, self.option("known_ko"))
        command = self.add_command("kegg_anno", cmd, ignore_error=True).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行kegg注释脚本完成")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("运行kegg注释脚本出错", code = "33701503")
        self.option('kegg_table', self.output_dir + '/kegg_table.xls')
