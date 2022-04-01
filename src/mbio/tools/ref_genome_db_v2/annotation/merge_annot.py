# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last update liubinxu 20180522
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


class MergeAnnotAgent(Agent):
    """
    将已知（参考基因组）序列和新序列的注释结果合一起
    """
    def __init__(self, parent):
        super(MergeAnnotAgent, self).__init__(parent)
        options = [
            {"name": "database", "type": "string", "default": "go,cog,kegg,pfam,nr,swissprot"},
            {"name": "annot_class", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "annot_class_ref", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "annot_db", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "is_assemble", "type": "bool", "default": True},
        ]
        self.add_option(options)
        self._memory_increase_step = 50
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
            raise OptionError("至少选择一种注释库", code = "33701901")
        if self.option("is_assemble") == True:
            if not self.option("annot_class").is_set:
                raise OptionError("需要指定新转录本注释路径", code = "33701902")
            if not self.option("annot_db").is_set:
                raise OptionError("需要指定新转录本数据库比对路径", code = "33701903")

        for db in self.database:
            if db not in ["go", "cog", "kegg", "pfam", "stat", "swissprot", "nr"]:
                raise OptionError("需要合并的注释文件不在支持范围内", code = "33701904")
            # if db == "go" and not self.option("gos_dir"):
            #     raise OptionError("缺少go注释的输入文件目录")
            # if db == "cog" and not self.option("cog_table_dir"):
            #     raise OptionError("缺少cog注释table的输入文件目录")
            # if db == "kegg":
            #     if not self.option("kegg_table_dir") and not self.option("pathway_table_dir"):
            #        raise OptionError("缺少kegg注释table和pathway_table的输入文件目录")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        self._memory = "50G"

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
        self.python = "/program/Python/bin/python"
        self.merge_scripts = self.config.PACKAGE_DIR + "/ref_rna_v2/merge.py"
        self.goAnnot = self.config.PACKAGE_DIR + "/ref_rna_v2/goAnnot.py"
        self.goSplit = self.config.PACKAGE_DIR + "/ref_rna_v2/goSplit.py"
        self.map_path = self.config.PACKAGE_DIR + "/ref_rna_v2/map4.r"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        if os.path.exists("/usr/bin/convert"):
            self.image_magick = "/usr/bin/convert"
        else:
            self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.merge_kegg_pathway = self.config.PACKAGE_DIR + "/ref_rna_v2/merge_kegg_pathway.py"
        self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"


        self.dir_dict = {
            "all_annot_tab_both" : "anno_stat/all_anno_detail.xls",
            "all_stat_tab_both" : "anno_stat/all_annotation_statistics.xls",
            "all_tran2gene_list_both"  :"tran2gene.txt",
            "cog_list_tab_tran" : "cog/cog.xls",
            "cog_summary_tab_gene" : "anno_stat/cog_stat/gene_cog_summary.xls",
            "cog_summary_tab_tran" : "cog/cog_summary.xls",
            "cog_venn_list_gene" : "anno_stat/venn/gene_cog_venn.txt",
            "cog_venn_list_tran" : "anno_stat/venn/cog_venn.txt",
            "go_lev2_stat_gene" : "anno_stat/go_stat/gene_go12level_statistics.xls",
            "go_lev2_stat_tran" : "go/go12level_statistics.xls",
            "go_lev3_stat_gene" : "anno_stat/go_stat/gene_go123level_statistics.xls",
            "go_lev3_stat_tran" : "go/go123level_statistics.xls",
            "go_lev4_stat_gene" : "anno_stat/go_stat/gene_go1234level_statistics.xls",
            "go_lev4_stat_tran" : "go/go1234level_statistics.xls",
            "go_list_tab_gene" : "anno_stat/go_stat/gene_gos.list",
            "go_list_tab_tran" : "go/query_gos.list",
            "go_venn_list_gene" : "anno_stat/venn/gene_go_venn.txt",
            "go_venn_list_tran" : "anno_stat/venn/go_venn.txt",
            "kegg_blast_tab_gene" : "anno_stat/blast/gene_kegg.xls",
            "kegg_blast_xml_gene" : "anno_stat/blast/gene_kegg.xml",
            "kegg_gene_tab_gene" : "anno_stat/kegg_stat/gene_kegg_table.xls",
            "kegg_gene_tab_tran" : "kegg/kegg_table.xls",
            "kegg_layer_tab_gene" : "anno_stat/kegg_stat/gene_kegg_layer.xls",
            "kegg_layer_tab_tran" : "kegg/kegg_layer.xls",
            "kegg_pathway_dir_gene" : "anno_stat/kegg_stat/gene_pathway",
            "kegg_pathway_dir_tran" : "kegg/pathways",
            "kegg_pathway_tab_gene" : "anno_stat/kegg_stat/gene_pathway_table.xls",
            "kegg_pathway_tab_tran" : "kegg/pathway_table.xls",
            "kegg_venn_list_gene" : "anno_stat/venn/gene_kegg_venn.txt",
            "kegg_venn_list_tran" : "anno_stat/venn/kegg_venn.txt",
            "nr_blast_tab_gene" : "anno_stat/blast/gene_nr.xls",
            "nr_blast_tab_tran" : "anno_stat/blast/nr.xls",
            "nr_blast_xml_gene" : "anno_stat/blast/gene_nr.xml", # no need
            "nr_venn_list_gene" : "anno_stat/venn/gene_nr_venn.txt",
            "nr_venn_list_tran" : "anno_stat/venn/nr_venn.txt",
            "pfam_domain_tab_gene" : "anno_stat/pfam_stat/gene_pfam_domain",
            "pfam_domain_tab_tran" : "pfam_domain",
            "pfam_venn_list_gene" : "anno_stat/venn/gene_pfam_venn.txt",
            "pfam_venn_list_tran" : "anno_stat/venn/pfam_venn.txt",
            "swissprot_blast_tab_gene" : "anno_stat/blast/gene_swissprot.xls",
            "swissprot_blast_tab_tran" : "anno_stat/blast/swissprot.xls",
            "swissprot_blast_xml_gene" : "anno_stat/blast/gene_swissprot.xml",
            "swissprot_evalue_tab_gene" : "anno_stat/blast_swissprot_statistics/gene_swissprot_evalue.xls",
            "swissprot_evalue_tab_tran" : "anno_stat/blast_swissprot_statistics/swissprot_evalue.xls",
            "swissprot_similar_tab_gene" : "anno_stat/blast_swissprot_statistics/gene_swissprot_similar.xls",
            "swissprot_similar_tab_tran" : "anno_stat/blast_swissprot_statistics/swissprot_similar.xls",
            "swissprot_venn_list_gene" : "anno_stat/venn/gene_swissprot_venn.txt",
            "swissprot_venn_list_tran" : "anno_stat/venn/swissprot_venn.txt"
        }

    def get_db_annot(self):
        '''
        获取数据库注释路径信息
        ; 返回参考注释路径
        '''
        f = open(self.json_path, "r")
        json_dict = json.loads(f.read())
        if json_dict.has_key(self.option("species_name")):
            return os.path.join(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/",
                                json_dict[self.option("species_name")]['anno_path_v2'])
        else:
            self.logger.info("数据库中没有该物种{}".format(self.option("species_name")))

    def creat_output(self):
        '''
        创建目录
        '''
        if os.path.exists(self.output_dir + "/newannot_mapdb"):
            shutil.rmtree(self.output_dir + "/newannot_mapdb")
        if self.option("is_assemble") == True:
            CopyFile().linkdir(self.option("annot_db").prop['path'], self.output_dir + "/newannot_mapdb")
        for annot_class in ["newannot_class", "refannot_class", "allannot_class"]:
            if os.path.exists(self.output_dir + "/" + annot_class):
                shutil.rmtree(self.output_dir + "/" + annot_class)
            os.makedirs(self.output_dir + "/" + annot_class )
            for db in self.option("database").split(","):
                os.makedirs(self.output_dir + "/" + annot_class + "/" + db)

    def link_result(self, result_dir="", out_dir=""):
        '''
        移动文件到结果目录中
        '''
        tail_dict = {"tab":"xls", "list":"txt", "xml":"xml", "stat":"stat.xls", "dir":"dir"}
        for file_key,file_dir in self.dir_dict.items():
            if os.path.exists(result_dir + "/" + file_dir):
                db,analysis,file_type,seq_type = file_key.split("_")
                if db == "all":
                    file_to_dir = os.path.join(self.output_dir, out_dir,
                                               "{}_{}.{}".format(db, analysis, tail_dict[file_type]))
                else:
                    file_to_dir = os.path.join(self.output_dir, out_dir, db,
                                               "{}_{}_{}.{}".format(db, analysis, seq_type, tail_dict[file_type]))
                if file_type == "dir":
                    file_to_dir = os.path.join(self.output_dir, out_dir, db,
                                               "{}_{}_{}_{}".format(db, analysis, seq_type, tail_dict[file_type]))
                    CopyFile().linkdir(result_dir + "/" + file_dir, file_to_dir)
                else:
                    CopyFile().linkfile(result_dir + "/" + file_dir, file_to_dir)
            else:
                self.logger.info("注释结果中没有该文件{}".format(file_dir))

    def run_merge_all(self, result_ref, result_new):
        '''
        合并注释结果
        '''
        tail_dict = {"tab":"xls", "list":"txt", "xml":"xml", "stat":"stat.xls", "dir":"dir"}
        for file_key,file_dir in self.dir_dict.items():
            if os.path.exists(result_ref + "/" + file_dir) and os.path.exists(result_new + "/" + file_dir):
                db,analysis,file_type,seq_type = file_key.split("_")
                if db == "all":
                    file_to_dir = os.path.join(self.output_dir, "allannot_class",
                                               "{}_{}.{}".format(db, analysis, tail_dict[file_type]))
                else:
                    file_to_dir = os.path.join(self.output_dir, "allannot_class", db,
                                               "{}_{}_{}.{}".format(db, analysis, seq_type, tail_dict[file_type]))
                if file_type == "dir":
                    file_to_dir = os.path.join(self.output_dir, "allannot_class", db,
                                               "{}_{}_{}_{}".format(db, analysis, seq_type, tail_dict[file_type]))
                self.merge2file(file_key, result_ref + "/" + file_dir, result_new + "/" + file_dir, file_to_dir)
            else:
                self.logger.info("注释结果中没有该文件{}".format(file_dir))

    def merge2file(self, file_key, file1, file2, file_to):
        '''
        合并两个文件到结果中
        '''
        if file_key in ["all_tran2gene_list_both",
                        "cog_venn_list_gene", "cog_venn_list_tran",
                        "go_list_tab_gene", "go_list_tab_tran", "go_venn_list_gene", "go_venn_list_tran",
                        "kegg_venn_list_gene", "kegg_venn_list_tran",
                        "nr_venn_list_gene", "nr_venn_list_tran",
                        "pfam_venn_list_gene", "pfam_venn_list_tran",
                        "swissprot_venn_list_gene", "swissprot_venn_list_tran"]:
            Merge().merge2file(file1=file1, file2=file2, fileto=file_to, withhead=False)
        elif file_key in ["all_annot_tab_both", "cog_list_tab_tran",
                          "nr_blast_tab_gene", "nr_blast_tab_tran",
                          "kegg_gene_tab_gene", "kegg_gene_tab_tran",
                          "pfam_domain_tab_gene", "pfam_domain_tab_tran",
                          "swissprot_blast_tab_gene", "swissprot_blast_tab_tran"]:
            Merge().merge2file(file1=file1, file2=file2, fileto=file_to, withhead=True)
        else:
            pass

        if file_key == "all_stat_tab_both":
            # 合并注释统计表
            pd1=pd.read_table(file1)
            pd2=pd.read_table(file2)
            pd3 = pd.DataFrame()
            pd3['type'] = pd1['type']
            pd3['transcripts'] = pd1['transcripts'] + pd2['transcripts']
            pd3['genes'] = pd1['genes'] + pd2['genes']
            pd3['transcripts_percent'] = pd3['transcripts']/list(pd3['transcripts'])[-1]
            pd3['genes_percent'] = pd3['genes']/list(pd3['genes'])[-1]
            pd3.to_csv(file_to, index=False, sep='\t', header=True)
        if file_key == "go_list_tab_gene":
            self.run_go_anno(file_to, "gene")
            os.system("mv go12level_statistics.xls {}/allannot_class/go/go_lev2_gene.stat.xls".format(self.output_dir))
            os.system("mv go123level_statistics.xls {}/allannot_class/go/go_lev3_gene.stat.xls".format(self.output_dir))
            os.system("mv go1234level_statistics.xls {}/allannot_class/go/go_lev4_gene.stat.xls".format(self.output_dir))
        if file_key == "go_list_tab_tran":
            self.run_go_anno(file_to, "tran")
            os.system("mv go12level_statistics.xls {}/allannot_class/go/go_lev2_tran.stat.xls".format(self.output_dir))
            os.system("mv go123level_statistics.xls {}/allannot_class/go/go_lev3_tran.stat.xls".format(self.output_dir))
            os.system("mv go1234level_statistics.xls {}/allannot_class/go/go_lev4_tran.stat.xls".format(self.output_dir))

    def merge_kegg(self):
        self.run_kegg_anno(self.output_dir + "/allannot_class/kegg/kegg_pathway_gene.xls",
                           self.option("annot_class_ref").prop['path'] + "/" + self.dir_dict["kegg_pathway_tab_gene"],
                           self.option("annot_class").prop['path'] + "/" + self.dir_dict["kegg_pathway_tab_gene"],
                           "gene",
                           self.output_dir + "/allannot_class/kegg/kegg_pathway_gene_dir",
                           self.output_dir + "/allannot_class/kegg/kegg_gene_gene.xls",
                          )
        self.run_kegg_anno(self.output_dir + "/allannot_class/kegg/kegg_pathway_tran.xls",
                           self.option("annot_class_ref").prop['path'] + "/" + self.dir_dict["kegg_pathway_tab_tran"],
                           self.option("annot_class").prop['path'] + "/" + self.dir_dict["kegg_pathway_tab_tran"],
                           "tran",
                           self.output_dir + "/allannot_class/kegg/kegg_pathway_tran_dir",
                           self.output_dir + "/allannot_class/kegg/kegg_gene_tran.xls",
                           )


    def merge_go(self):
        pass

    def run_kegg_anno(self, kegg_table, r_level_path, n_level_path, seq_type, path_out, merged_gene_kegg):
        # 合并kegg注释
        cmd = "{} {} {} {} {} {} {} {} {} {} {}".format(
            self.python,
            self.merge_kegg_pathway,
            self.map_path,
            self.r_path,
            self.image_magick,
            r_level_path,
            n_level_path,
            kegg_table,
            path_out,
            merged_gene_kegg,
            self.html_path
        )
        self.logger.info("开始画图")
        self.logger.info(cmd)
        cmd1_obj = self.add_command("merge_png" + seq_type, cmd, ignore_error=True).run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("画图完成")
        elif cmd1_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("画图出错", code = "33701905")



    def run_go_anno(self, go_list, seq_type):
        cmd1 = "{} {} {} {} {} {}".format(self.python, self.goAnnot, go_list, "localhost", self.b2g_user, self.b2g_password)
        self.logger.info("运行goAnnot.py")
        self.logger.info(cmd1)
        cmd1_obj = self.add_command("cmd1" + seq_type, cmd1, ignore_error=True).run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("运行goAnnot.py完成")
        elif cmd1_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd1_obj.rerun()
            self.wait(cmd1_obj)
            if cmd1_obj.return_code == 0:
                self.logger.info("运行goAnnot.py完成")
            else:
                self.set_error("运行goAnnot.py出错", code = "33701906")

        cmd2 = "{} {} {}".format(self.python, self.goSplit, self.work_dir + '/go_detail.xls')
        self.logger.info("运行goSplit.py")
        self.logger.info(cmd2)
        cmd2_obj = self.add_command("cmd2" + seq_type, cmd2, ignore_error=True).run()
        self.wait(cmd2_obj)
        if cmd2_obj.return_code == 0:
            self.logger.info("运行goSplit.py完成")
        elif cmd2_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd2_obj.rerun()
            self.wait(cmd2_obj)
            if cmd2_obj.return_code == 0:
                self.logger.info("运行goSplit.py完成")
            else:
                self.set_error("运行goSplit.py出错", code = "33701907")


    def merge(self, dirs, merge_file, type):
        cmd = "{} {} {} {}".format(self.python, self.merge_scripts, dirs, merge_file)
        cmd3_obj = self.add_command("cmd_merge_{}".format(type), cmd, ignore_error=True).run()
        self.wait(cmd3_obj)
        if cmd3_obj.return_code == 0:
            self.logger.info("文件合并完成")
        elif cmd3_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd3_obj.rerun()
            self.wait(cmd3_obj)
            if cmd3_obj.return_code == 0:
                self.logger.info("文件合并完成")
            else:
                self.set_error("文件合并未完成", code = "33701908")


    def run(self):
        super(MergeAnnotTool, self).run()
        self.creat_output()
        # ref_annot = self.get_db_annot()
        if self.option("annot_class_ref").is_set:
            self.ref_path = self.option("annot_class_ref").prop['path']
        if self.option("annot_class_ref").is_set:
            self.link_result(self.ref_path, "refannot_class")
        if self.option("is_assemble") == True:
            self.new_path = self.option("annot_class").prop['path']
            self.link_result(self.new_path, "newannot_class")
        if self.option("annot_class_ref").is_set and self.option("annot_class").is_set:
            self.run_merge_all(self.ref_path, self.new_path)
            self.merge_kegg()
        self.end()
