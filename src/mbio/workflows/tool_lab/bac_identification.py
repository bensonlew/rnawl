# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
from biocluster.workflow import Workflow
import datetime
import types
import shutil
from bson.objectid import ObjectId
import gevent
from mbio.packages.metagbin.common_function import link_dir,link_file
import HTMLParser
from mbio.packages.tool_lab.common import down_seq_files


class BacIdentificationWorkflow(Workflow):
    """
    细菌菌种鉴定小工具
    """
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        super(BacIdentificationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'},
            {'name': "input_format", 'type': 'string'},
            {'name': 'input_dir', 'type': 'infile', 'format': 'toolapps.fasta_dir'},
            {'name': 'input_file', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'gff_dir', 'type': 'infile', 'format': 'tool_lab.gff_dir'},
            {'name': 'gff_file', 'type': 'infile', 'format': 'gene_structure.gff3'},
            {'name': "ref_database", 'type': 'string', "default": "GTDB"},  ##数据库GTDB or custom
            {'name': "ref_sample_list", 'type': 'string'},  ###"sample1;sample2;sample3"
            {'name': 'project_data', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.genome_pre = self.add_module('toolapps.genomes_predict')
        self.gff_genome_pre = self.add_module('toolapps.gff_genome_predict')
        self.ncbi_download = self.add_module('toolapps.ncbi_download')
        self.genome_ani = self.add_module('toolapps.genomes_taxon')
        if self.option('project_data'):
            print(self._sheet)
            self.project_data = eval(HTMLParser.HTMLParser().unescape(self._sheet.option('project_data')))
            for i in self.project_data['specimens']:
                self.samples[i['id']] = i['name']
            assemble_dir = os.path.join(self._sheet.work_dir, "assemble_dir")
            if os.path.exists(assemble_dir):
                shutil.rmtree(assemble_dir)
            (self.assemble_dir, self.analysis_type) = down_seq_files(self.project_data['my_type'],
                                                                     self.project_data['db_version'], assemble_dir,
                                                                     self._sheet.option("project_task_id"),
                                                                     self.samples)
        
    def run(self):
        if self.option("project_data"):
            self.genome =self.download_file()
            if self.option("ref_database") in ["GTDB"]:
                self.genome_pre.on("end", self.run_genome_taxon)
                self.genome_ani.on("end", self.set_output)
                self.run_genomes()
            elif self.option("ref_database") in ["custom"]:
                self.on_rely([self.genome_pre, self.ncbi_download], self.run_genome_taxon)
                self.genome_ani.on("end", self.set_output)
                self.run_genomes()
                self.run_custom()
        else:
            self.get_dir()
            if len(os.listdir(self.work_dir + "/gff")) == 0 and len(
                    os.listdir(self.work_dir + "/genome")) > 0 and self.option("ref_database") in ["GTDB"]:
                self.genome_pre.on("end", self.run_genome_taxon)
                self.genome_ani.on("end", self.set_output)
                self.run_genomes()
            elif len(os.listdir(self.work_dir + "/gff")) > 0 and len(
                    os.listdir(self.work_dir + "/genome")) > 0 and self.option("ref_database") in ["GTDB"]:
                self.gff_genome_pre.on("end", self.run_genome_taxon)
                self.genome_ani.on("end", self.set_output)
                self.run_gff_genomes()
            elif len(os.listdir(self.work_dir + "/gff")) == 0 and len(
                    os.listdir(self.work_dir + "/genome")) > 0 and self.option("ref_database") in ["custom"]:
                self.on_rely([self.genome_pre, self.ncbi_download], self.run_genome_taxon)
                self.genome_ani.on("end", self.set_output)
                self.run_genomes()
                self.run_custom()
            elif len(os.listdir(self.work_dir + "/gff")) > 0 and len(
                    os.listdir(self.work_dir + "/genome")) > 0 and self.option("ref_database") in ["custom"]:
                self.on_rely([self.gff_genome_pre, self.ncbi_download], self.run_genome_taxon)
                self.genome_ani.on("end", self.set_output)
                self.run_gff_genomes()
                self.run_custom()
        super(BacIdentificationWorkflow, self).run()


    def get_dir(self):
        if os.path.exists(self.work_dir+"/genome"):
            shutil.rmtree(self.work_dir+"/genome")
        os.mkdir(self.work_dir+"/genome")
        if os.path.exists(self.work_dir+"/gff"):
            shutil.rmtree(self.work_dir+"/gff")
        os.mkdir(self.work_dir+"/gff")
        if self.option("input_file").is_set and not self.option("gff_file").is_set:
            sample = ".".join(os.path.basename(self.option("input_file").prop['path']).split(".")[0:-1])
            os.link(self.option("input_file").prop['path'], self.work_dir+"/genome/"+sample+".fasta")
        elif self.option("input_file").is_set and self.option("gff_file").is_set:
            sample = ".".join(os.path.basename(self.option("input_file").prop['path']).split(".")[0:-1])
            os.link(self.option("input_file").prop['path'], self.work_dir + "/genome/" + sample + ".fasta")
            os.link(self.option("gff_file").prop['path'], self.work_dir + "/gff/" + sample + ".gff")
            with open(self.work_dir+"/mapping_file.xls", "w") as g:
                g.write("{}\t{}\t{}\n".format(sample, sample + ".fasta", sample + ".gff"))
        if self.option("input_dir").is_set and not self.option("gff_dir").is_set:
            with open(self.option("input_dir").prop['path']+"/list.txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    lin = line.strip().split("\t")
                    os.link(lin[0], self.work_dir + "/genome/" + lin[1] + ".fasta")
        elif self.option("input_dir").is_set and self.option("gff_dir").is_set:
            genome_dict ={}
            with open(self.option("input_dir").prop['path'] + "/list.txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    lin = line.strip().split("\t")
                    os.link(lin[0], self.work_dir + "/genome/" + lin[1] + ".fasta")
                    genome_dict[lin[1]] = lin[1] + ".fasta"
            gff_dict ={}
            with open(self.option("gff_dir").prop['path']+ "/list.txt", "r") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    lin = line.strip().split("\t")
                    self.logger.info(self.option("gff_dir").prop['path'])
                    os.link(self.option("gff_dir").prop['path']+"/"+lin[0], self.work_dir + "/gff/" + lin[1] + ".gff")
                    gff_dict[lin[1]] = lin[1] + ".gff"
            with open(self.work_dir+"/mapping_file.xls", "w") as g:
                for key in genome_dict.keys():
                    if key in gff_dict.keys():
                        g.write("{}\t{}\t{}\n".format(key, genome_dict[key], gff_dict[key]))

    def download_file(self):
        """
        download file from s3
        :return:
        """
        path = ''
        if self.analysis_type in ['complete']:
            assemble_dir2 = os.path.join(self.work_dir, "assemble_dir2")
            if os.path.exists(assemble_dir2):
                shutil.rmtree(assemble_dir2)
            os.mkdir(assemble_dir2)
            for sample in self.samples.keys():
                sample_path = os.path.join(assemble_dir2, self.samples[sample] + ".fasta")
                assemble_dir3 = os.path.join(self.assemble_dir, self.samples[sample])
                dir_list = os.listdir(assemble_dir3)
                for file2 in dir_list:
                    os.system("cat {} >> {}".format(os.path.join(assemble_dir3, file2), sample_path))
            path = assemble_dir2
        elif self.analysis_type in ['uncomplete']:
            for sample in self.samples.keys():
                sample_path = os.path.join(assemble_dir2, self.samples[sample] + ".fasta")
                assemble_dir3 = os.path.join(self.assemble_dir, self.samples[sample])
        return path

    def run_genomes(self):
        """
        只运行只有基因组序列的情况(单个文件或多个文件)
        :return:
        """
        genome = ''
        if self.option("project_data"):
            genome = self.genome
        else:
            genome = self.work_dir + "/genome"
        opts = {
            "genomes": genome
            }
        self.genome_pre.set_options(opts)
        self.genome_pre.run()

    def run_gff_genomes(self):
        """
        只运行只有基因组序列和gff的情况(单个文件或多个文件)
        :return:
        """
        opts = {
            "genomes": self.work_dir + "/genome",
            "gffs": self.work_dir + "/gff",
            "map_info": self.work_dir+"/mapping_file.xls",
        }
        self.gff_genome_pre.set_options(opts)
        self.gff_genome_pre.run()

    def run_custom(self):
        """
        只运行自定义数据库
        :return:
        """
        opts = {
            "sample_list": self.option("ref_sample_list"),
        }
        self.ncbi_download.set_options(opts)
        self.ncbi_download.run()

    def run_genome_taxon(self):
        """
        基因组物种注释分类
        :return:
        """
        opts = ''
        house_dir = ''
        s16_dir = ''
        genome_fa = ''
        if self.option("project_data"):
            genome_fa =self.genome
            s16_dir = self.genome_pre.output_dir + "/s16"
            house_dir = self.genome_pre.output_dir + "/house"
        else:
            genome_fa = self.work_dir + "/genome"
            if len(os.listdir(self.work_dir + "/gff")) == 0 and len(os.listdir(self.work_dir + "/genome")) > 0:
                s16_dir = self.genome_pre.output_dir + "/s16"
                house_dir = self.genome_pre.output_dir + "/house"
            elif len(os.listdir(self.work_dir + "/gff")) > 0 and len(os.listdir(self.work_dir + "/genome")) > 0:
                s16_dir = self.gff_genome_pre.output_dir + "/s16"
                house_dir = self.gff_genome_pre.output_dir + "/house"
        if self.option("ref_database") in ["custom"]:
            custom_fa = self.ncbi_download.option("s16")
            cus_table = self.ncbi_download.option("cus_table")
            opts = {
                "type": self.option("ref_database"),
                "16s": s16_dir,
                "house": house_dir,
                "genomes": genome_fa,
                "ref_genomes":self.ncbi_download.option("genomes"),
                "ref_s16": self.ncbi_download.option("s16"),
                "custom_fa": custom_fa,
                "cus_table": cus_table
            }
        elif self.option("ref_database") in ["GTDB"]:
            opts = {
                "type": self.option("ref_database"),
                "16s": s16_dir,
                "house": house_dir,
                "genomes": genome_fa,
            }
        self.genome_ani.set_options(opts)
        self.genome_ani.run()

    def set_output(self):
        """
        结果文件
        :return:
        """
        if self.option("project_data"):
            if self.option("ref_database") in ["GTDB"]:
                link_file(self.genome_pre.option("stat").prop['path'], self.work_dir + "/all.stat.xls")
                link_file(self.genome_ani.output_dir + "/all.taxon.xls", self.work_dir + "/all.taxon.xls")
                if len(os.listdir(self.genome_ani.output_dir + "/anno")) > 0:
                    link_dir(self.genome_ani.output_dir + "/anno", self.output_dir + "/anno")
                else:
                    os.mkdir(self.output_dir + "/anno")
                if len(os.listdir(self.genome_ani.output_dir + "/tree")) > 0:
                    link_dir(self.genome_ani.output_dir + "/tree", self.output_dir + "/tree")
                else:
                    if os.path.exists(self.output_dir + "/tree"):
                        shutil.rmtree(self.output_dir + "/tree")
                    os.mkdir(self.output_dir + "/tree")
                if len(os.listdir(self.genome_ani.output_dir + "/16s_blast")) > 0:
                    link_dir(self.genome_ani.output_dir + "/16s_blast", self.output_dir + "/16s_blast")
            elif self.option("ref_database") in ["custom"]:
                if os.path.exists(self.work_dir + "/all.stat.xls"):
                    os.remove(self.work_dir + "/all.stat.xls")
                os.system("cat {} {} > {}".format(self.genome_pre.option("stat").prop['path'],
                                                  self.ncbi_download.option("database_stat").prop['path'],
                                                  self.work_dir + "/all.stat.xls"))
                link_file(self.genome_ani.output_dir + "/all.taxon.xls", self.work_dir + "/all.taxon.xls")
                if len(os.listdir(self.genome_ani.output_dir + "/anno")) > 0:
                    link_dir(self.genome_ani.output_dir + "/anno", self.output_dir + "/anno")
                else:
                    os.mkdir(self.output_dir + "/anno")
                if len(os.listdir(self.genome_ani.output_dir + "/tree")) > 0:
                    link_dir(self.genome_ani.output_dir + "/tree", self.output_dir + "/tree")
                else:
                    if os.path.exists(self.output_dir + "/tree"):
                        shutil.rmtree(self.output_dir + "/tree")
                    os.mkdir(self.output_dir + "/tree")
                if len(os.listdir(self.genome_ani.output_dir + "/16s_blast")) > 0:
                    link_dir(self.genome_ani.output_dir + "/16s_blast", self.output_dir + "/16s_blast")
        else:
            if len(os.listdir(self.work_dir + "/gff")) == 0 and len(
                    os.listdir(self.work_dir + "/genome")) > 0 and self.option("ref_database") in ["GTDB"]:
                link_file(self.genome_pre.option("stat").prop['path'], self.work_dir + "/all.stat.xls")
                link_file(self.genome_ani.output_dir + "/all.taxon.xls", self.work_dir + "/all.taxon.xls")
                if len(os.listdir(self.genome_ani.output_dir + "/anno")) > 0:
                    link_dir(self.genome_ani.output_dir + "/anno", self.output_dir + "/anno")
                else:
                    os.mkdir(self.output_dir + "/anno")
                if len(os.listdir(self.genome_ani.output_dir + "/tree")) > 0:
                    link_dir(self.genome_ani.output_dir + "/tree", self.output_dir + "/tree")
                else:
                    if os.path.exists(self.output_dir + "/tree"):
                        shutil.rmtree(self.output_dir + "/tree")
                    os.mkdir(self.output_dir + "/tree")
                if len(os.listdir(self.genome_ani.output_dir + "/16s_blast")) > 0:
                    link_dir(self.genome_ani.output_dir + "/16s_blast", self.output_dir + "/16s_blast")
            elif len(os.listdir(self.work_dir + "/gff")) > 0 and len(
                    os.listdir(self.work_dir + "/genome")) > 0 and self.option("ref_database") in ["GTDB"]:
                link_file(self.gff_genome_pre.option("stat").prop['path'], self.work_dir + "/all.stat.xls")
                link_file(self.genome_ani.output_dir + "/all.taxon.xls", self.work_dir + "/all.taxon.xls")
                if len(os.listdir(self.genome_ani.output_dir + "/anno")) > 0:
                    link_dir(self.genome_ani.output_dir + "/anno", self.output_dir + "/anno")
                else:
                    os.mkdir(self.output_dir + "/anno")
                if len(os.listdir(self.genome_ani.output_dir + "/tree")) > 0:
                    link_dir(self.genome_ani.output_dir + "/tree", self.output_dir + "/tree")
                else:
                    if os.path.exists(self.output_dir + "/tree"):
                        shutil.rmtree(self.output_dir + "/tree")
                    os.mkdir(self.output_dir + "/tree")
                if len(os.listdir(self.genome_ani.output_dir + "/16s_blast")) > 0:
                    link_dir(self.genome_ani.output_dir + "/16s_blast", self.output_dir + "/16s_blast")
            elif len(os.listdir(self.work_dir + "/gff")) == 0 and len(
                    os.listdir(self.work_dir + "/genome")) > 0 and self.option("ref_database") in ["custom"]:
                if os.path.exists(self.work_dir + "/all.stat.xls"):
                    os.remove(self.work_dir + "/all.stat.xls")
                os.system("cat {} {} > {}".format(self.genome_pre.option("stat").prop['path'],
                                                  self.ncbi_download.option("database_stat").prop['path'],
                                                  self.work_dir + "/all.stat.xls"))
                link_file(self.genome_ani.output_dir + "/all.taxon.xls", self.work_dir + "/all.taxon.xls")
                if len(os.listdir(self.genome_ani.output_dir + "/anno")) > 0:
                    link_dir(self.genome_ani.output_dir + "/anno", self.output_dir + "/anno")
                else:
                    os.mkdir(self.output_dir + "/anno")
                if len(os.listdir(self.genome_ani.output_dir + "/tree")) > 0:
                    link_dir(self.genome_ani.output_dir + "/tree", self.output_dir + "/tree")
                else:
                    if os.path.exists(self.output_dir + "/tree"):
                        shutil.rmtree(self.output_dir + "/tree")
                    os.mkdir(self.output_dir + "/tree")
                if len(os.listdir(self.genome_ani.output_dir + "/16s_blast")) > 0:
                    link_dir(self.genome_ani.output_dir + "/16s_blast", self.output_dir + "/16s_blast")
            elif len(os.listdir(self.work_dir + "/gff")) > 0 and len(
                    os.listdir(self.work_dir + "/genome")) > 0 and self.option("ref_database") in ["custom"]:
                if os.path.exists(self.work_dir + "/all.stat.xls"):
                    os.remove(self.work_dir + "/all.stat.xls")
                os.system("cat {} {} > {}".format(self.gff_genome_pre.option("stat").prop['path'],
                                                  self.ncbi_download.option("database_stat").prop['path'],
                                                  self.work_dir + "/all.stat.xls"))
                link_file(self.genome_ani.output_dir + "/all.taxon.xls", self.work_dir + "/all.taxon.xls")
                if len(os.listdir(self.genome_ani.output_dir + "/anno")) > 0:
                    link_dir(self.genome_ani.output_dir + "/anno", self.output_dir + "/anno")
                else:
                    if os.path.exists(self.output_dir + "/anno"):
                        shutil.rmtree(self.output_dir + "/anno")
                    os.mkdir(self.output_dir + "/anno")
                if len(os.listdir(self.genome_ani.output_dir + "/tree")) > 0:
                    link_dir(self.genome_ani.output_dir + "/tree", self.output_dir + "/tree")
                else:
                    if os.path.exists(self.output_dir + "/tree"):
                        shutil.rmtree(self.output_dir + "/tree")
                    os.mkdir(self.output_dir + "/tree")
                if len(os.listdir(self.genome_ani.output_dir + "/16s_blast")) > 0:
                    link_dir(self.genome_ani.output_dir + "/16s_blast", self.output_dir + "/16s_blast")
        self.get_files(self.work_dir, self.output_dir)
        self.set_db()
        self.end()

    def set_db(self):
        self.logger.info("导表开始")
        api_bac = self.api.api("tool_lab.bac_identif")
        stat_file =self.work_dir + "/all.stat.xls"
        taxon_file = self.work_dir + "/all.taxon.xls"
        dir = self.output_dir+"/tree"
        dir2 = self.output_dir+"/anno"
        main_id = ObjectId(self.option("main_id"))
        if self.option("ref_database") in ["GTDB"]:
            api_bac.add_detail(main_id, stat_file, taxon_file, dir, dir2)
        elif self.option("ref_database") in ["custom"]:
            api_bac.add_detail(main_id, stat_file, taxon_file, dir, dir2, undefined="custom")
        self.logger.info("导表结束")
        if len(os.listdir(self.output_dir + "/tree")) == 0:
            shutil.rmtree(self.output_dir + "/tree")
        if len(os.listdir(self.output_dir + "/anno")) == 0:
            shutil.rmtree(self.output_dir + "/anno")

    def get_files(self, dir , dir2):
        if os.path.exists(dir2+"/all.taxon.xls"):
            os.remove(dir2+"/all.taxon.xls")
        if os.path.exists(dir2+"/all.stat.xls"):
            os.remove(dir2+"/all.stat.xls")
        if self.option("ref_database") in ["custom"]:
            with open(dir+"/all.taxon.xls", "r") as f,open(dir2+"/all.taxon.xls", "w") as g:
                lines = f.readlines()
                g.write("Sample\tMost Similar taxonomy（NCBI)\n")
                for line in lines:
                    g.write(line)
            with open(dir+"/all.stat.xls", "r") as f,open(dir2+"/all.stat.xls", "w") as g:
                lines = f.readlines()
                g.write("Sample\t16S_num\thousekeeping_num\tType\n")
                for line in lines:
                    g.write(line)
        elif self.option("ref_database") in ["GTDB"]:
            with open(dir + "/all.taxon.xls", "r") as f, open(dir2 + "/all.taxon.xls", "w") as g:
                lines = f.readlines()
                g.write("Sample\tMost Similar taxonomy（GTDB)\n")
                for line in lines:
                    g.write(line)
            with open(dir + "/all.stat.xls", "r") as f, open(dir2 + "/all.stat.xls", "w") as g:
                lines = f.readlines()
                g.write("Sample\t16S_num\thousekeeping_num\n")
                for line in lines:
                    g.write(line)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "菌种鉴定分析目录", 0],
            ["./all.taxon.xls", "xls", "所有样品菌种鉴定汇总表", 0],
            ["./all.stat.xls", "txt", "所有样品的信息统计表", 0],
            ["./anno/", "", "所有样品注释的目录", 0],
            ["./16s_blast/", "", "所有样品16s blast结果目录", 0],
            ["./tree/", "", "所有样品的进化树目录", 0],
        ])
        result_dir.add_regexp_rules([
            [r'.house_keeping.nwk', 'nwk', '看家基因进化树', 0],
            [r'.16s.nwk', 'nwk', '16s进化树', 0],
            [r'.16s.blast.xls', 'xls', '16s比对结果', 0],
            [r'.anno.xls', 'xls', '单个样品的注释详情表', 0],
        ])
        super(BacIdentificationWorkflow, self).end()