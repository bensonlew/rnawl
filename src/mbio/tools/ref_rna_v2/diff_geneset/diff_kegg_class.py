# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from collections import defaultdict
import gridfs
import re
from bson import ObjectId
import time
from mbio.packages.itraq_and_tmt.kegg_regulate import KeggRegulate
# from mbio.packages.ref_rna_v2.kegg_regulate import KeggRegulate
from biocluster.config import Config
import subprocess
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
from mbio.packages.rna.annot_config import AnnotConfig
import shutil
import unittest
from mbio.packages.ref_rna_v2.functions import tryforgood

class DiffKeggClassAgent(Agent):
    """
    Kegg分类统计分析，主要用于基因集的重运行步骤
    version v1.0.1
    author: shijin
    last_modify: 2017.8.16
    """
    def __init__(self, parent):
        super(DiffKeggClassAgent, self).__init__(parent)
        options = [
            {"name": "geneset_kegg", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "task_id", "type": "string"},
            {"name": "kegg_table", "type": "infile", "format": "medical_transcriptome.kegg_table"},
            {"name": "kegg_table2", "type": "string", "default": ""},
            {"name": "regulate", "type": "string", "default": ""},
            {"name": "geneset_id", "type": "string"},
            {"name": "level", "type": "string"},
            {"name": "background_links", "type": "string", "default": ""},  # 底图的地址信息，add_info
            # {"name": "type", "type": "string", "default": "origin"},  # 取最新的注释表还是原来的注释表
            {"name": "source", "type": "string", "default": ""},
            {'name': 'kegg_version', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.step.add_steps("kegg_regulate")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kegg_regulate.start()
        self.step.update()

    def stepfinish(self):
        self.step.kegg_regulate.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('{} = {}'.format("kegg_table", self.option("kegg_table").prop["path"]))
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 11
        self._memory = '10G'

    def end(self):
        super(DiffKeggClassAgent, self).end()


class DiffKeggClassTool(Tool):
    def __init__(self, config):
        super(DiffKeggClassTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python = '/program/Python/bin/'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map5.r"

        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))

        self.db_path = self.kegg_files_dict["html"]
        self.db_path_old = self.config.SOFTWARE_DIR + "/database/KEGG/xml/"
        # if self.option('kegg_version') in ["201909", "202003"]:
        #     self.db_path = os.path.join(self.config.SOFTWARE_DIR,
        #                                 'database/Annotation/other{}/kegg{}/html/'.format(self.option('kegg_version')[0:4],
        #                                 self.option('kegg_version')))

        # else:
        #     self.db_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"

        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.species_abr = 'map'
        self.gene2set = dict()
        self.geneset_ko = list()
        self.geneset_list = list()
        self.genome_path = self.config.SOFTWARE_DIR + "/database/KEGG/genome2.xls"
        self.map_dict = {}
        # self.geneset_id_number = len(self.option('geneset_id').split(','))
        # self.colors = ['#DC143C', '#008000'] if self.option('source') == 'diff_exp' else ["#0000cd", "#ff8000"]
        # self.colors = ['#FF0000', '#00FF00']
        if self.option("regulate") in ["all","up"]:
            self.colors = ['#FF0000', '#0000FF']
        else:
            self.colors = [ '#0000FF','#FF0000']
        self.genome_path = self.kegg_files_dict["genome2.xls"]

    def run(self):
        """
        运行
        :return:
        """
        super(DiffKeggClassTool, self).run()
        # self.get_kegg_pics()
        self.get_kegg_png()
        self.get_dicts()
        self.get_kegg_stat_tmp()
        if self.option("background_links"):
            self.get_background_info()
        self.generate_kegg_stat()
        pathways = self.output_dir + '/pathways'
        if not os.path.exists(pathways):
            os.mkdir(pathways)
        self.generate_ko_txt_dir2()
        self.generate_new_pics()
        self.packing_press()
        self.end()

    def get_genome_abr(self, kegg_species):
        with open(self.genome_path, 'r' ) as f:
            lines = f.readlines()
            for line in lines:
                genome_id = re.sub("gn:", "", line.split("\t")[0])
                genome_abr = line.split("\t")[1].split(',')[0]
                genome = line.split("\t")[1].split(';')[-1].strip()
                if kegg_species == genome_abr or kegg_species == genome:
                    self.species_abr = genome_abr
                    break

    def get_kegg_pics(self):
        """该函数被get_kegg_png替换"""
        self.mongo_db = Config().get_mongo_client(mtype="medical_transcriptome")[Config().get_mongo_dbname("medical_transcriptome")]
        fs = gridfs.GridFS(self.mongo_db)
        annotation_collection = self.mongo_db["sg_annotation_kegg"]
        geneset_collection = self.mongo_db["sg_geneset"]
        anno_type = self.option("level")
        main_id = annotation_collection.find_one({"task_id": self.option('task_id')})["main_id"]
        kegg_level_collection = self.mongo_db["sg_annotation_kegg_level"]
        results = kegg_level_collection.find({"kegg_id":main_id, "anno_type":anno_type})
        anno_path = self.work_dir + "/png"
        if not os.path.exists(anno_path):
            os.mkdir(anno_path)
        self.logger.info("开始导出kegg png格式图片")
        start = time.time()
        for result in results:
            graph_id = result["graph_png_id"]
            pathway_id = result["pathway_id"]
            with open(anno_path + "/" + pathway_id + ".png", "w") as fw:
                content = fs.get(graph_id).read()
                fw.write(content)
        end = time.time()
        self.logger.info("png格式图片生成完成，耗时{}".format(end - start))

    @tryforgood
    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            if os.path.isfile(to_path) or os.path.islink("raw_png"):
                os.remove(to_path)
            else:
                try:
                    shutil.rmtree(to_path)
                except:
                    os.system('rm -r {}'.format(to_path.rstrip("/")))
        if os.path.exists(path):
            os.system('ln -s {} {}'.format(path.rstrip("/"), to_path.rstrip("/")))
        else:
            try:
                transfer = MultiFileTransfer()
                transfer.add_download(path, to_path)
                transfer.perform()
            except:
                return False
        return to_path

    def get_kegg_png(self):

        # 工作流跑的时候还需要判定一下物种
        import pandas as pd
        kegg_df = pd.read_csv(self.option('kegg_table').prop['path'], sep='\t').fillna('')
        paths = sum([p.split(';') for p in kegg_df['Paths'] if p], [])
        specie = paths[0]
        self.species_abr = filter(str.isalpha, specie)
        if self.species_abr.lower() == 'ko':
            self.species_abr = 'map'

        self.mongo_db = Config().get_mongo_client(mtype="medical_transcriptome")[Config().get_mongo_dbname("medical_transcriptome")]
        kegg_collection = self.mongo_db["sg_annotation_stat"]
        geneset_collection = self.mongo_db["sg_geneset"]

        # task_id = result["task_id"]
        geneset_type = self.option("level")
        kegg_annot_output = kegg_collection.find_one({"task_id": self.option('task_id'), 'status': "end"})['result_dir']

        class_dir = 'allannot_class'
        if kegg_collection.find_one({"task_id": self.option('task_id'), 'status': "end"})['has_new'] == False:
            class_dir = 'refannot_class'
        if geneset_type == 'G':
            kegg_pathways = os.path.join(kegg_annot_output, class_dir, 'kegg', 'kegg_pathway_gene_dir')
        else:
            kegg_pathways = os.path.join(kegg_annot_output, class_dir, 'kegg', 'kegg_pathway_tran_dir')
        anno_path = self.work_dir + "/png/"

        if not kegg_pathways.endswith("/"):
            kegg_pathways += "/"
        anno_path = self.work_dir + "/raw_png/"
        self.logger.info("lalalalalalal{}".format(kegg_pathways))
        try:
            self.download_s3_file(kegg_pathways, anno_path)
        except:
            self.set_error('file can not find %s', variables=kegg_pathways, code="33706404")
        if os.path.exists("png"):
            shutil.rmtree("png")
        if not os.path.exists("raw_png"):
            self.set_error('raw_png下载失败')
        shutil.copytree("raw_png", "png")
        # os.system('ln -s {} {}'.format(kegg_pathways, anno_path))

    def get_dicts(self):
        gene_kegg_gene = dict()
        if self.species_abr == 'map':
            ko_genes, path_ko = self.option('kegg_table').get_pathway_koid( kegg_version=self.option('kegg_version'))
        else:
            #区分比对到整个ko与单个物种的ko时有区别
            ko_genes, path_ko, gene_kegg_gene = self.option('kegg_table').get_pathway_koid2(kegg_version=self.option('kegg_version'))
        geneset_ko = defaultdict(set)
        regulate_gene = {}
        gene2set = dict()
        geneset2_ko = list()
        with open(self.option("geneset_kegg"), "r") as f:
            for line in f:
                line = line.strip().split("\t")
                regulate_gene[line[0]] = line[1].split(",")
                geneset_ko[line[0]] = []
                self.geneset_list.append(line[0])
                for key in ko_genes.keys():
                    for gene in regulate_gene[line[0]]:
                        if gene in ko_genes[key]:
                            geneset_ko[line[0]].append(key)
                geneset2_ko.append(set(geneset_ko[line[0]]))
                for gene in regulate_gene[line[0]]:
                    if gene2set.has_key(gene):
                        gene2set[gene].append(line[0])
                    else:
                        gene2set[gene] = [line[0]]

                    if self.species_abr == 'map':
                        pass
                    else:
                        if gene_kegg_gene.has_key(gene):
                            for kegg_genes in gene_kegg_gene[gene]:
                                for kegg_gene in kegg_genes.split(";"):
                                    geneset_ko[line[0]].append(kegg_gene)
                        else:
                            pass
        for i in gene2set.keys():
            self.gene2set[i] = ",".join(gene2set[i])
        if self.gene2set.has_key(""):
            del self.gene2set[""]
        self.logger.info("path is {}".format(";".join(path_ko.keys())))
        self.logger.info("gene2set is {}".format(self.gene2set))
        self.ko_genes= ko_genes  # ko与基因的对应关系
        self.category = geneset_ko  # 基因集与ko的对应关系
        self.path_ko = path_ko  # path小ko与基因大Ko的对应关系
        self.geneset_gene = regulate_gene  # 基因集与基因的对应关系
        self.geneset_ko = geneset2_ko
        self.gene_kegg_gene = gene_kegg_gene


        ## 根据基因集过滤通路
        # geneset_map_list = list()
        # with open(self.option('kegg_table').prop['path'], 'r') as r:
        #     r.readline()
        #     for line in r:
        #         line = line.strip('\n').split('\t')
        #         if line[0] in gene2set:
        #             paths = [re.sub('path:', '', i) for i in line[-1].split(';')]
        #             geneset_map_list.extend(paths)

        # self.path_ko = {k: v for k,v in path_ko.items() if k in geneset_map_list}
        # self.path_ko = path_ko  # path小ko与基因大Ko的对应关系



        self.geneset_gene = regulate_gene  # 基因集与基因的对应关系
        self.gene_kegg_gene = gene_kegg_gene #基因与kegg基因的对应关系

    def get_kegg_stat_tmp(self):
        self.logger.info("生成kegg_stat_tmp.xls文件")
        if self.species_abr == 'map':
            KeggRegulate().get_regulate_table(ko_gene=self.ko_genes, path_ko=self.path_ko,
                                          regulate_gene=self.geneset_gene,
                                          output= self.work_dir + '/kegg_stat_tmp.xls',
                                          gene_list = self.geneset_list)
        else:
            KeggRegulate().get_regulate_table2(ko_gene=self.ko_genes, path_ko=self.path_ko,
                                               regulate_gene=self.geneset_gene,
                                               gene_kegg_gene=self.gene_kegg_gene,
                                               output=self.work_dir + '/kegg_stat_tmp.xls',
                                               gene_list =self.geneset_list)
        self.logger.info("生成kegg_stat_tmp.xls文件完毕")
    def generate_ko_txt_dir2(self):
        kegg_regulate=self.output_dir + '/kegg_stat.xls'
        out_dir=self.output_dir
        if not os.path.exists(out_dir + "/ko"):
            os.mkdir(out_dir + "/ko")
        f = open(kegg_regulate, "r")
        f.readline()
        for line in f:
            tmp = line.strip().split("\t")
            path = tmp[0]
            gene_num1 = tmp[2]
            gene1_list = []
            gene2_list = []

            # 提取注释到该通路的基因集列表
            if len(self.geneset_gene) == 1:
                gene1_list = [x.split("(")[0] for x in tmp[3].split(");")]
            else:
                if gene_num1 and not gene_num1 == "0" and not gene_num1.startswith("http"):
                    gene1_list = [x.split("(")[0] for x in tmp[3].split(");")]
                gene_num2 = tmp[4]
                if gene_num2 and not gene_num2.startswith("http") and not gene_num2 == "0":
                    gene2_list = [x.split("(")[0] for x in tmp[5].split(");")]
            # 提取颜色信息
            ko_list_unrepeat = set(tmp[1].split(";"))
            color_dict = {}
            for ko in ko_list_unrepeat:
                ko_genes_set = set(self.ko_genes[ko])
                color_dict[ko] = []
                self.logger.info("ko{}\tset1{}\tset2{}".format(ko_genes_set, gene1_list, gene2_list))
                if len(self.geneset_gene) == 1:
                    if  ko_genes_set.intersection(gene1_list):
                        color_dict[ko].append("#0000cd")  # 蓝色
                elif len(self.geneset_gene) ==2:
                    if ko_genes_set.intersection(gene1_list) and ko_genes_set.intersection(gene2_list):
                        color_dict[ko].append("#0000cd,#ff0000")  # 蓝色,大红
                    elif ko_genes_set.intersection(gene1_list):
                        color_dict[ko].append("#0000cd")  # 蓝色
                    elif ko_genes_set.intersection(gene2_list):
                        color_dict[ko].append("#ff0000")  # 大红
                else:
                    pass
            with open(out_dir + "/ko/" + path, "w") as fw:
                fw.write("#KO\tbg\tfg\n")
                for key in color_dict.keys():
                    if len(color_dict[key]) != 0:
                        str_ = key + "\tNA\t" + ",".join(color_dict[key]) + "\n"
                    else:
                        str_ = key + "\tNA\t" + "NA" + "\n"
                    fw.write(str_)
        self.logger.info("ko文件生成完毕")
        return



    def generate_ko_txt_dir(self):
        # catgory=self.geneset_gene
        kegg_regulate=self.output_dir + '/kegg_stat.xls'
        out_dir=self.output_dir
        if not os.path.exists(out_dir + "/ko"):
            os.mkdir(out_dir + "/ko")
        f = open(kegg_regulate, "r")
        f.readline()
        for line in f:
            tmp = line.strip().split("\t")
            path = tmp[0]
            gene_num1 = tmp[2]
            gene1_list = []
            gene2_list = []
            if len(self.geneset_gene) == 1:
                if gene_num1 and not gene_num1 == "0" and not gene_num1.startswith("http"):
                    gene1_list_tmp = [x.split("(")[1] if not ")" in x else x.strip(")").split("(")[1] for x in tmp[3].split(");")]
                    for gene in gene1_list_tmp:
                        if gene.find(";") != 1:
                            gene1_list.append(gene)
                        else:
                            gene1_list.extend(gene.split(";"))
            else:
                if gene_num1 and not gene_num1 == "0" and not gene_num1.startswith("http"):
                    gene1_list_tmp = [x.split("(")[1] if not ")" in x else x.strip(")").split("(")[1] for x in tmp[3].split(");")]
                    for gene in gene1_list_tmp:
                        if gene.find(";") != 1:
                            gene1_list.append(gene)
                        else:
                            gene1_list.extend(gene.split(";"))
                gene_num2 = tmp[4]
                if gene_num2 and not gene_num2.startswith("http") and not gene_num2 == "0":
                    gene2_list_tmp = [x.split("(")[1] if not ")" in x else x.strip(")").split("(")[1] for x in tmp[5].split(");")]
                    for gene in gene2_list_tmp:
                        if gene.find(";") != 1:
                            gene2_list.append(gene)
                        else:
                            gene2_list.extend(gene.split(";"))
            gene_list = []
            gene_list.extend(gene1_list)
            gene_list.extend(gene2_list)
            gene_list_unrepeat = list(set(gene_list))
            color_dict = {}
            for gene in gene_list_unrepeat:
                color_dict[gene] = []
                if len(self.geneset_gene) == 1:
                    if gene in gene1_list:
                        color_dict[gene].append(self.colors[0])  # 蓝色
                elif len(self.geneset_gene) ==2:
                    if gene in gene1_list and gene in gene2_list:
                        color_dict[gene].append(','.join(self.colors))  # 蓝色,大红
                    elif gene in gene1_list:
                        color_dict[gene].append(self.colors[0])  # 蓝色
                    elif gene in gene2_list:
                        color_dict[gene].append(self.colors[1])  # 大红
                else:
                    pass
            with open(out_dir + "/ko/" + path, "w") as fw:
                fw.write("#KO\tbg\tfg\n")
                for key in color_dict.keys():
                    if len(color_dict[key]) != 0:
                        str_ = key + "\tNA\t" + ",".join(color_dict[key]) + "\n"
                    else:
                        str_ = key + "\tNA\t" + "NA" + "\n"
                    fw.write(str_)
        self.logger.info("ko文件生成完毕")
        return


    def generate_new_pics(self):
        path_ko=self.path_ko
        kos_path=self.output_dir + "/ko"
        out_dir=self.output_dir
        path_list = path_ko.keys()
        cmd1_list = []
        cmd2_list = []
        for path in path_list:
            ko_path = kos_path + "/" + path
            if os.path.exists(ko_path):
                ko_old = path.replace(self.species_abr, "ko")
                ko = path.replace(self.species_abr, "map")
                if ko.startswith("ko"):
                    ko.replace("ko", "map")
                ko_xml = self.db_path + "/" + ko + ".kgml"
                if self.option("kegg_version") <= "2020":
                    if os.path.exists(self.db_path_old + ko_old + ".xml"):
                        ko_xml = self.db_path_old + ko_old + ".xml"
                    else:
                        ko_xml = self.db_path + ko + ".kgml"

                cmd = "{} {} {} {} {} {} {}".format(self.r_path, self.map_path, path,
                                                 ko_path, out_dir + "/pathways/" + path + ".png",
                                                 ko_xml,
                                                 self.work_dir + "/png/" + path + ".png")
                map_html = KeggHtml(version=self.option('kegg_version'))
                if self.option("regulate") in ["all", "up"]:
                    map_html.color_fg = ['#FF0000', '#0000FF']
                else:
                    map_html.color_fg = ['#0000FF', '#FF0000']
                # try:
                #     map_html.run_gene_set(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html", self.gene2set)
                # except:
                #     self.logger.info("注释结果中没有kegghtml {}".format(path))
                if os.path.exists(self.work_dir + "/png/" + path + ".html.mark"):
                    map_html.run_gene_set_mark(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/png/" + path + ".html.mark", self.geneset_ko)
                elif os.path.exists(self.work_dir + "/png/" + path + ".html"):
                    map_html.run_gene_set_mark_from_html(path + '.png', self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/png/" + path + ".html", self.geneset_ko)
                elif os.path.exists(self.work_dir + "/kegg/" + path + ".html.mark"):
                    map_html.run_gene_set_mark(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/kegg/" + path + ".html.mark", self.geneset_ko)
                else:
                    self.logger.info("{}注释结果中没有html 或 mark文件".format(path))
                db_png_path = self.work_dir + "/png/" + path + ".png"
                if os.path.exists(out_dir + "/pathways/" + path + ".png"):
                    os.remove(out_dir + "/pathways/" + path + ".png")
                if os.path.exists(db_png_path):
                    os.link(db_png_path, out_dir + "/pathways/" + path + ".png")
                cmd1_list.append(cmd)
                ## 去掉pnd转成pdf的步骤，运行太耗时
                # pdf_path = out_dir + "/pathways/" + path + ".pdf"
                # cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                #       + out_dir + "/pathways/" + path + ".png" + ' ' + pdf_path
                # cmd2_list.append(cmd)
            else:
                if os.path.exists(self.work_dir + "/png/"):
                    db_png_path =  self.work_dir + "/png/" + path + ".png"
                else:
                    db_png_path = self.work_dir + "/kegg/" + path + ".png"
                if os.path.exists(out_dir + "/pathways/" + path + ".png"):
                    os.remove(out_dir + "/pathways/" + path + ".png")
                if os.path.exists(db_png_path):
                    os.link(db_png_path, out_dir + "/pathways/" + path + ".png")
                else:
                    print(db_png_path)
                ## 去掉pnd转成pdf的步骤，运行太耗时
                # cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                #       + db_png_path + ' ' + out_dir + "/pathways/" + path + ".pdf"
                # cmd1_list.append(cmd)
        self.logger.info("开始生成新kegg图片")
        with open(self.work_dir + "/cmd1.list", "w") as fw:
            for i in range(len(cmd1_list)):
                fw.write(cmd1_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd1.list")
        cmd1_obj = self.add_command("cmd1", cmd, ignore_error=True).run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("cmd1 list执行成功")
        with open(self.work_dir + "/cmd2.list", "w") as fw:
            for i in range(len(cmd2_list)):
                fw.write(cmd2_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd2.list")
        cmd2_obj = self.add_command("cmd2", cmd, ignore_error=True).run()
        self.wait(cmd2_obj)
        if cmd2_obj.return_code == 0:
            self.logger.info("cmd2 list执行成功")
        self.logger.info("kegg图片生成完毕")


    def get_background_info(self):
        background_links = self.option("background_links")
        if not os.path.exists(background_links):
            self.logger.info("不存在背景信息文件，程序退出")
            return
        with open(background_links) as r:
            r.readline()
            for line in r:
                tmp = line.strip().split("\t")
                pathway = tmp[0]
                links = tmp[1].split("?")[1].split("/")
                links.pop(0)  # 去除掉map
                dict = {}
                for link in links:
                    ko = link.split("%09")[0]  # K06101
                    color = link.split("%09")[1]  # tomato
                    dict[ko] = color
                self.map_dict[pathway] = dict  # self.map_dict["map05340"] = {"K10887": "yellow"}

    def packing_press(self):
        if os.path.exists(os.path.join(self.output_dir,"pathways.tar.gz")):
            os.remove(os.path.join(self.output_dir,"pathways.tar.gz"))
        cmd = "tar -zcvf {} {} ".format(os.path.join(self.output_dir,"pathways.tar.gz"), os.path.join(self.output_dir,"pathways"))
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('kegg图片打包压缩完成')
        except subprocess.CalledProcessError:
            self.logger.info('kegg图片打包压缩失败')
            self.set_error("kegg图片打包压缩失败失败")

    def generate_kegg_stat(self):
        with open(self.work_dir + '/kegg_stat_tmp.xls', "r") as file, open(self.output_dir + '/kegg_stat.xls', "w") as fw:
            line = file.readline()
            fw.write(line)
            for line in file:
                tmp = line.strip().split("\t")
                links = tmp[-1]
                pathway = tmp[0]
                links_ko = links.split("?")[1].split("/")
                links_ko.pop(0)  # 去除掉map
                lnk = links.split("?")[0] + "?map=" + pathway + "&multi_query="
                ko_tmp = [x.split("%09")[0] for x in links_ko]
                ko_tmp2 = []
                for x in ko_tmp:
                    ko_tmp2.extend(x.split(';'))
                if pathway in self.map_dict:  # 含有背景色
                    for ko in self.map_dict[pathway].keys():  # 对背景色中的所有项进行循环
                        # self.logger.info(ko)
                        if ko == "":
                            continue
                        if ko in ko_tmp2:  # 基因ko显著富集
                            # self.logger.info("{} {}".format(ko, ",".join(ko_tmp2)))
                            if self.get_color(ko):
                                lnk += ko + "+{},{}%0d%0a".format(self.map_dict[pathway][ko], self.get_color(ko))  # 有背景色,前景色
                            else:
                                lnk += ko + "+{}%0d".format(self.map_dict[pathway][ko])  # 只有背景色
                        else:
                            lnk += ko + "+{}%0d".format(self.map_dict[pathway][ko])  # 只有背景色
                else:  # 只标边框
                    for ko in ko_tmp2:
                        if ko == "":
                            continue
                        if self.get_color(ko):
                            lnk += ko + "+{},{}%0d%0a".format("white", self.get_color(ko))  # 只有背景色
                        else:
                            pass
                tmp[-1] = lnk
                str_ = "\t".join(tmp) + "\n"
                fw.write(str_)

    def get_color(self, ko):
        """
        根据基因的ko号，获取颜色
        :param ko: 基因的ko号
        :return:
        """
        if len(self.category) == 1:
            if ko in self.category[self.category.keys()[0]]:
                return "red"
                # return "blue"
            else:
                return False
        elif len(self.category) == 2:
            # if self.option('source') == 'diff_exp':
            # lst = list(self.category.keys())  # 基因集列表
            lst = self.geneset_list
            # lst.sort()
            if ko in self.category[lst[0]] and ko in self.category[lst[1]]:
                return "pink"
            elif ko in self.category[lst[0]]:
                return "red"
            elif ko in self.category[lst[1]]:
                return "blue"
            else:
                return False

            '''
            else:
                pass

            # lst = list(self.category.keys())  # 基因集列表
            lst = self.geneset_list
            # lst.sort()
            if ko in self.category[lst[0]] and ko in self.category[lst[1]]:
                return "pink"
            elif ko in self.category[lst[0]]:
                return "blue"
            elif ko in self.category[lst[1]]:
                return "red"
            else:
                return False
            '''

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "test_keggcalssref" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "medical_transcriptome.diff_geneset.diff_kegg_class",
            "instant": False,
            "options": dict(
                task_id = "emvhefm6875rraot9k7c0jkv70",
                geneset_kegg="/mnt/ilustre/users/sanger-dev/workspace/20201228/DiffGenesetPipline_emvhefm6875rraot9k7c0jkv70_8343_2519/DiffGenesetPrepare/DiffGenesetPrepare1/output/kegg_class/multi_geneset_list",
                #group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
                level="G",
                kegg_table="/mnt/ilustre/users/sanger-dev/workspace/20201228/DiffGenesetPipline_emvhefm6875rraot9k7c0jkv70_8343_2519/DiffGenesetPrepare/AnnotPrepare/output/gene_kegg_table.xls",
                kegg_table2="/mnt/ilustre/users/sanger-dev/workspace/20201228/DiffGenesetPipline_emvhefm6875rraot9k7c0jkv70_8343_2519/DiffGenesetPrepare/AnnotPrepare/output/gene_kegg_level_table.xls",
                background_links="/mnt/ilustre/users/sanger-dev/workspace/20201228/DiffGenesetPipline_emvhefm6875rraot9k7c0jkv70_8343_2519/DiffGenesetPrepare/AnnotPrepare/output/add_info.txt",
                source="diff_exp",
                # annot_result = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test20201110/annot",
                kegg_version = "202007"

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
