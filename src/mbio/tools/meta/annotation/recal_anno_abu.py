# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180417
# last modified : guhaidong

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
from  mainapp.models.mongo.metagenomic import  Metagenomic
from biocluster.config import Config


#from mbio.packages.annotation.mg_annotation.kegg_pathway_img import KeggPathwayImg

class RecalAnnoAbuAgent(Agent):
    """
    宏基因组注释交互分析功能筛选
    """

    def __init__(self, parent):
        super(RecalAnnoAbuAgent, self).__init__(parent)
        options = [
            {"name": "gene_anno_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "database", "type": "string"},  # nr、cog、kegg、cazy、ardb、card、vfdb
            {"name": "origin_anno", "type": "infile", "format": "sequence.profile_table"},  # kegg时使用
            {"name": "xml_file", "type": "infile", "format": "sequence.profile_table"},  # kegg时使用
            {"name": "anno_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"},
            {"name": "vfdb_type", "type": "string", "default": "core"},  # vfdb时使用
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # vfdb时使用
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老任务
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # modified by qingchen.zhang@ 20191118

    def check_options(self):
        if not self.option("gene_anno_table").is_set:
            raise OptionError("必须设置基因注释文件", code="32700801")
        if not self.option("gene_profile").is_set:
            raise OptionError("必须设置基因丰度文件", code="32700802")
        '''
        if not self.option("database") in ["nr", "cog", "kegg", "cazy", "ardb", "card", "vfdb"]:
            raise OptionError("请输入正确的数据库名称", code="32700803")
        '''
        if not self.option("origin_anno").is_set and self.option("database") == "kegg":
            raise OptionError("输入kegg时必须输入origin_anno文件!", code="32700804")
        # if not self.option("xml_file").is_set and self.option("database") == "kegg":
            # raise OptionError("输入kegg时必须输入xml_file文件!", code="32700805")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'  # 改回 by guhaidong @ 20180427
        # memory = 5 + 10 * self._rerun_time  # 每次重运行增加5G内存 by guhaidong @ 20180417
        # self._memory = "%sG" % memory

    def end(self):
        self.option("anno_result_dir", self.output_dir)
        super(RecalAnnoAbuAgent, self).end()


class RecalAnnoAbuTool(Tool):
    def __init__(self, config):
        super(RecalAnnoAbuTool, self).__init__(config)
        self.python_path = "/program/Python/bin/python"
        self.perl_path = self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.nr_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/nr_anno_abudance.pl'
        self.cog_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/eggNOG_anno_abundance.pl'
        self.kegg_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/kegg_anno_abudance.pl'
        self.cazy_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/cazy_anno_abu.pl'


        self.vfdb_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/vfdb_anno_abu_new.pl'
        self.vfdb_des = self.config.SOFTWARE_DIR + "/database/CAZyDB/FamInfo.txt"
        self.probio_scr = self.config.PACKAGE_DIR + '/annotation/mg_annotation/probio_abu.pl'
        self.qs_scr = self.config.PACKAGE_DIR + '/annotation/qs_anno_abundance.pl'
        #self.go_scr = self.config.PACKAGE_DIR + '/annotation/qs_anno_abundance.pl'
        self.p450_scr = self.config.PACKAGE_DIR + '/annotation/mg_annotation/cyps_anno_abundance.pl'
        self.pfam_scr = self.config.PACKAGE_DIR + '/annotation/mg_annotation/pfam_anno_abundance.pl'

        self.metagenomic = Metagenomic() ## add by qingchen.zhang ## 兼容不同版本的数据库
        self.metagenomic._config = Config()
        self.version_value = ''
        self.version_value = self.metagenomic.find_version_from_task(self.option("task_id"))
        if self.version_value in [""]:
            self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
            self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/kegg_pathway_img.py'
            self.card_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/card_anno_abudance.pl'
            self.ardb_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/ardb_anno_abudance_old.pl' ##
            # ardb区分新老版本，原因是新老任务表头不一样，可能会导致后面的筛选出现问题
        else: ## 94.2
            self.card_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/card_anno_abudance_v3.0.9.pl'
            self.html_path = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/html/"
            self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/kegg_pathway_img_v94.py'
            self.ardb_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/ardb_anno_abudance.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(RecalAnnoAbuTool, self).run()
        self.run_calculate_abu()
        if self.option('database') == "kegg":
            self.run_kegg_img()
        self.set_output()
        self.end()

    def run_calculate_abu(self):
        self.logger.info("start recalculate abundance")
        geneprofile = self.option('gene_profile').prop['path']
        gene_anno = self.option('gene_anno_table').prop['path']
        if self.option("database") == "nr":
            cmd = '{} {} -q {} -p {} -o {}'. \
                format(self.perl_path, self.nr_script, gene_anno, geneprofile, self.output_dir)
        elif self.option("database") == "cog":
            cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.cog_script, gene_anno,
                                                   geneprofile, self.output_dir)
        elif self.option("database") == "kegg":
            self.logger.info(gene_anno)
            origin_anno = self.option('origin_anno').prop['path']
            paths = origin_anno.split("/")
            anno_dir = "/".join(paths[0:len(paths) - 1])
            self.logger.info(anno_dir)
            e_file = os.path.join(anno_dir, "kegg_enzyme_profile.xls")
            m_file = os.path.join(anno_dir, "kegg_module_profile.xls")
            p_file = os.path.join(anno_dir, "kegg_pathway_profile.xls")
            cmd = "{} {} -q {} -p {} -e {} -m {} -path {} -o {}".format(self.perl_path, self.kegg_script,
                                                                        gene_anno, geneprofile, e_file, m_file,
                                                                        p_file, self.output_dir)
        elif self.option("database") == "cazy":
            cmd = "{} {} -q {} -p {} -des {} -o {}".format(self.perl_path, self.cazy_script, gene_anno,
                                                           geneprofile, self.vfdb_des, self.output_dir)
        elif self.option("database") == "ardb":
            cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.ardb_script, gene_anno,
                                                   geneprofile, self.output_dir)
        elif self.option("database") == "card":
            cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.card_script, gene_anno,
                                                   geneprofile, self.output_dir)
        elif self.option("database") == "vfdb":
            type = self.option("vfdb_type")
            if self.option("group").is_set:
                group = self.option("group").prop["path"]
                cmd = "{} {} -c {} -p {} -type {} -group {} -o {}".format(self.perl_path, self.vfdb_script, gene_anno,
                                                                          geneprofile, type, group, self.output_dir)
            else:
                cmd = "{} {} -c {} -p {} -type {} -o {}".format(self.perl_path, self.vfdb_script, gene_anno,
                                                                geneprofile, type, self.output_dir)
        elif self.option("database") == "probio":
            cmd = "{} {} -i {} -p {} -o {}".format(self.perl_path, self.probio_scr, gene_anno,
                                                   geneprofile, self.output_dir)
        elif self.option("database") == "qs":
            cmd = "{} {} -i {} -geneprofile {} -o {}".format(self.perl_path, self.qs_scr, gene_anno,
                                                   geneprofile, self.output_dir)
        elif self.option("database") == "p450":
            cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.p450_scr, gene_anno,
                                                   geneprofile, self.output_dir)
        elif self.option("database") == "pfam":
            cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.pfam_scr, gene_anno,
                                                   geneprofile, self.output_dir)
        command = self.add_command("recalculate_abundance", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("recalculate abundance succeed")
        elif command.return_code == -9:  # add memory limit by guhaidong @ 20180417
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("recalculate abundance failed", code="32700801")
            self.set_error("recalculate abundance failed", code="32700803")

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        try:
            self.option("anno_result_dir", self.output_dir)
            self.logger.info("设置输出结果目录成功")
        except Exception as e:
            self.set_error("输出结果文件目录异常——%s", variables=(e), code="32700804")

    def run_kegg_img(self):
        self.logger.info("start output kegg img")
        #run_img = KeggPathwayImg()
        # xml_file = self.option('xml_file').prop['path']
        pathway_file = self.output_dir + "/kegg_pathway_eachmap.xls"
        # cmd2 = "{} {} -i {} -o {} -p {} -ko {} -KO {}".format(self.python_path, script, xml_file, self.output_dir,
        #                                                       pathway_file, "Pathway", "KO_list")
        # cmd2 = "{} {} -i {} -o {} -p {} -ko {} -KO {} -png_file {} -html {}".format(self.python_path, self.script, xml_file, self.output_dir,pathway_file, "Pathway", "KO_list", "True", self.html_path)
        cmd2 = "{} {} -o {} -p {} -ko {} -KO {} -png_file {} -html {}".format(self.python_path, self.script,
                                                                                    self.output_dir,
                                                                                    pathway_file, "Pathway", "KO_list",
                                                                                    "True", self.html_path)
        self.logger.info(cmd2)
        command = self.add_command("output_kegg_pathway_img", cmd2, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("output_kegg_pathway_img succeed")
        else:
            self.set_error("output kegg pathway img failed", code="32700802")
            #run_img.run_img(xml_file, pathway_file, "#Pathway", "KO_list", self.output_dir)
