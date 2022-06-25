# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from collections import defaultdict
import gridfs
from bson import ObjectId
import time
import json
from mbio.packages.itraq_and_tmt.kegg_regulate import KeggRegulate
# from mbio.packages.rna.kegg_html import KeggHtml
from mbio.packages.protein_transcript.kegg_html import KeggHtml
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.config import Config
import subprocess
import tarfile
import shutil
from mbio.packages.rna.annot_config import AnnotConfig



class KeggClassAgent(Agent):
    """
    Kegg分类统计分析，主要用于基因集的重运行步骤
    version v1.0.1
    author: shijin
    last_modify: 2017.8.16
    """
    def __init__(self, parent):
        super(KeggClassAgent, self).__init__(parent)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "proteinset_kegg", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "kegg_table", "type": "infile", "format": "itraq_and_tmt.kegg_table"},
            {"name": "kegg_id", "type": "string"},
            {"name": "ko_txt", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': "2017"},
            {"name": "gene_kegg_table", "type": "string"},
            {"name": "background_links", "type": "string", "default": ""},  # 底图的地址信息，add_info
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
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 11
        self._memory = '10G'

    def end(self):
        super(KeggClassAgent, self).end()


class KeggClassTool(Tool):
    def __init__(self, config):
        super(KeggClassTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python = '/miniconda2/bin/'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map5.r"
        self.map_re = self.config.PACKAGE_DIR + "/protein_transcript/map4.r"
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))
        self.db_path_old = self.config.SOFTWARE_DIR + "/database/KEGG/xml/"
        self.db_path = self.kegg_files_dict["html"]
        # self.db_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
        self.geneset_ko = list()
        self.proteinset_list = list()
        self.map_dict = {}

        self.species_abr = 'map'
        self.gene2set = dict()
        self.genome_path = self.kegg_files_dict["genome2.xls"]
        # self.genome_path = self.config.SOFTWARE_DIR + "/database/KEGG/genome2.xls"
        self.k2g = dict()

        with open(self.option("gene_kegg_table")) as gr:
            _ = gr.readline()
            for line in gr:
                if line.strip():
                    line = line.strip().split('\t')
                    ks = line[1].split(";")
                    for k in ks:
                        if k not in self.k2g:
                            self.k2g[k] = set()
                        self.k2g[k].add(line[0])


    def run(self):
        """
        运行
        :return:
        """
        super(KeggClassTool, self).run()
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
        self.end()

    def get_genome_abr(self, kegg_species):
        with open (self.genome_path, 'r' ) as f:
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
        kegg_id = self.option("kegg_id")
        self.mongo_db = Config().get_mongo_client(mtype="itraq_tmt")[Config().get_mongo_dbname("itraq_tmt")]
        fs = gridfs.GridFS(self.mongo_db)
        annotation_collection = self.mongo_db["sg_annotation_kegg"]
        result = annotation_collection.find_one({"main_id": ObjectId(kegg_id)})
        if not result:
            self.set_error('kegg_annot with main_id: %s was not found', variables = (kegg_id), code = "32503501")

        kegg_level_collection = self.mongo_db["sg_annotation_kegg_level"]
        results = kegg_level_collection.find({"kegg_id":kegg_id})
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

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            if to_path.endswith('/'):
                os.system('ln -s {} {}'.format(path, to_path.rstrip('/')))
        else:
            try:
                transfer = MultiFileTransfer()
                transfer.add_download(path, to_path)
                transfer.perform()
            except:
                self.set_error('file can not find {}'.format(path))
        return to_path
    def get_kegg_png(self):
        kegg_id = self.option("kegg_id")
        self.mongo_db = Config().get_mongo_client(mtype="itraq_tmt")[Config().get_mongo_dbname("itraq_tmt")]
        annotation_collection = self.mongo_db["sg_annotation_stat"]
        kegg_collection = self.mongo_db["sg_annotation_kegg"]
        result = kegg_collection.find_one({"main_id": ObjectId(kegg_id)})
        if not result:
            self.set_error('kegg_annot with main_id: %s was not found', variables=(kegg_id), code="32503501")
        task_id = result["task_id"]
        anno_type = result["type"]
        kegg_annot_output = annotation_collection.find_one({"task_id": task_id, "type": anno_type, 'status': "end"})['result_dir']
        annot_params = annotation_collection.find_one({"task_id": task_id, "type": anno_type, 'status': "end"})['params']
        params_dict = json.loads(annot_params)
        if params_dict.has_key('kegg_species'):
            self.get_genome_abr(params_dict['kegg_species'])

        kegg_created_ts = annotation_collection.find_one({"task_id": task_id, "type": anno_type, 'status': "end"})['created_ts']
        if kegg_created_ts < "2021-07-28 18:10:49":
            kegg_pathways = os.path.join(kegg_annot_output, 'kegg', 'pathways')
            if not kegg_pathways.endswith("/"):
                kegg_pathways += "/"
            anno_path = self.work_dir + "/raw_png/"
            self.download_s3_file(kegg_pathways, anno_path)
            shutil.copytree("raw_png", "png")
            if os.path.exists(anno_path + "pathways.tar.gz"):
                tar = tarfile.open(anno_path + "pathways.tar.gz")
                tar.extractall(path=self.work_dir + "/png")
            elif os.path.exists(anno_path):
                pass
            else:
                self.set_error("无法获得kegg注释结果图片信息", code="32503504")
        else:
            # 旧的结果中，pathways是个目录，缺点是文件多传输慢。如今改为打包形式(pathways.tar)，需要下载包并解包。
            kegg_pathways = os.path.join(kegg_annot_output, 'kegg', 'pathways.tar')
            try:
                transfer = MultiFileTransfer()
                transfer.add_download(kegg_pathways, self.work_dir)
                transfer.perform()
            except:
                pass
            if not os.path.exists(os.path.join(self.work_dir, 'pathways.tar')):
                kegg_pathways = os.path.join(kegg_annot_output, '3_anno_KEGG', 'pathways.tar')
                try:
                    transfer = MultiFileTransfer()
                    transfer.add_download(kegg_pathways, self.work_dir)
                    transfer.perform()
                except:
                    self.set_error('无法获得kegg注释结果图片信息,file can not find %s', variables=(kegg_pathways), code="32503503")
            if os.path.exists(os.path.join(self.work_dir, 'pathways.tar')):
                pathways_tar = tarfile.open(os.path.join(self.work_dir, 'pathways.tar'))
                pathways_tar.extractall(path=self.work_dir)
                os.rename(os.path.join(self.work_dir, 'pathways'),os.path.join(self.work_dir, 'png'))

    def get_dicts(self):
        gene_kegg_gene = dict()
        if self.species_abr == 'map':
            ko_genes, path_ko = self.option('kegg_table').get_pathway_koid(kegg_version=self.option('kegg_version'))
            # ko_genes, path_ko = self.option('kegg_table').get_pathway_koid()
        else:
            ko_genes, path_ko, gene_kegg_gene = self.option('kegg_table').get_pathway_koid2(kegg_version=self.option('kegg_version'))
            # ko_genes, path_ko, gene_kegg_gene = self.option('kegg_table').get_pathway_koid2()
        proteinset_ko = defaultdict(set)
        regulate_gene = {}
        gene2set = dict()
        proteinset2_ko = list()
        with open(self.option("proteinset_kegg"), "r") as f:
            for line in f:
                line = line.strip().split("\t")
                regulate_gene[line[0]] = line[1].split(",")
                proteinset_ko[line[0]] = []
                self.proteinset_list.append(line[0])
                for key in ko_genes.keys():
                    for gene in regulate_gene[line[0]]:
                        if gene in ko_genes[key]:
                            proteinset_ko[line[0]].append(key)
                proteinset2_ko.append(set(proteinset_ko[line[0]]))
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
                                    proteinset_ko[line[0]].append(kegg_gene)
                        else:
                            pass
        for i in gene2set.keys():
            self.gene2set[i] = ",".join(gene2set[i])
        self.logger.info("path is {}".format(";".join(path_ko.keys())))

        self.category = proteinset_ko  # 基因集与ko的对应关系
        self.ko_genes= ko_genes  # ko与基因的对应关系
        self.path_ko = path_ko  # path小ko与基因大Ko的对应关系
        self.proteinset_gene = regulate_gene  # 基因集与基因的对应关系
        self.proteinset_ko = proteinset2_ko
        self.gene_kegg_gene = gene_kegg_gene #基因与kegg基因的对应关系

    def get_kegg_stat_tmp(self):
        self.logger.info("生成kegg_stat_tmp.xls文件")
        if self.species_abr == 'map':
            KeggRegulate().get_regulate_table(ko_gene=self.ko_genes, path_ko=self.path_ko,
                                              regulate_gene=self.proteinset_gene,
                                              output= self.work_dir + '/kegg_stat_tmp.xls')
        else:
            KeggRegulate().get_regulate_table2(ko_gene=self.ko_genes, path_ko=self.path_ko,
                                              regulate_gene=self.proteinset_gene,
                                              gene_kegg_gene=self.gene_kegg_gene,
                                              output= self.work_dir + '/kegg_stat_tmp.xls')

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
            if len(self.proteinset_gene) == 1:
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
                if len(self.proteinset_gene) == 1:
                    if  ko_genes_set.intersection(gene1_list):
                        color_dict[ko].append("#0000cd")  # 蓝色
                elif len(self.proteinset_gene) ==2:
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
        # catgory=self.proteinset_gene
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
            if len(self.proteinset_gene) == 1:
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
                if len(self.proteinset_gene) == 1:
                    if gene in gene1_list:
                        color_dict[gene].append("#0000cd")  # 蓝色
                elif len(self.proteinset_gene) ==2:
                    if gene in gene1_list and gene in gene2_list:
                        color_dict[gene].append("#0000cd,#ff0000")  # 蓝色,大红
                    elif gene in gene1_list:
                        color_dict[gene].append("#0000cd")  # 蓝色
                    elif gene in gene2_list:
                        color_dict[gene].append("#ff0000")  # 大红
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
        cmd0_list = []
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

                png_ = ko_xml.replace('ko', 'map').split('.')[0] + '.png'
                cmd0 = "{} {} {} {} {} {} {}".format(self.r_path, self.map_re, path,
                                                 self.option('ko_txt'), self.work_dir + "/png/" + path + ".png",
                                                 ko_xml,
                                                 png_)

                cmd = "{} {} {} {} {} {} {}".format(self.r_path, self.map_path, path,
                                                 ko_path, out_dir + "/pathways/" + path + ".png",
                                                 ko_xml,
                                                 self.work_dir + "/png/" + path + ".png")
                map_html = KeggHtml(version=self.option('kegg_version'))
                # try:
                #     map_html.run_gene_set(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html", self.gene2set)
                # except:
                #     self.logger.info("注释结果中没有kegghtml {}".format(path))
                if os.path.exists(self.work_dir + "/png/" + path + ".html.mark"):
                    map_html.run_gene_set_mark(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/png/" + path + ".html.mark", self.proteinset_ko, self.k2g)
                elif os.path.exists(self.work_dir + "/kegg/" + path + ".html.mark"):
                    map_html.run_gene_set_mark(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/kegg/" + path + ".html.mark", self.proteinset_ko, self.k2g)
                else:
                    self.logger.info("{}注释结果中没有html 或 mark文件".format(path))
                # print(self.proteinset_ko)
                db_png_path = self.work_dir + "/png/" + path + ".png"
                if os.path.exists(out_dir + "/pathways/" + path + ".png"):
                    os.remove(out_dir + "/pathways/" + path + ".png")
                if os.path.exists(db_png_path):
                    os.link(db_png_path, out_dir + "/pathways/" + path + ".png")
                cmd0_list.append(cmd0)
                cmd1_list.append(cmd)
                pdf_path = out_dir + "/pathways/" + path + ".pdf"
                cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                      + out_dir + "/pathways/" + path + ".png" + ' ' + pdf_path
                cmd2_list.append(cmd)
            else:
                if os.path.exists(self.work_dir + "/png/"):
                    db_png_path =  self.work_dir + "/png/" + path + ".png"
                else:
                    db_png_path = self.work_dir + "/kegg/" + path + ".png"
                if os.path.exists(out_dir + "/pathways/" + path + ".png"):
                    os.remove(out_dir + "/pathways/" + path + ".png")
                if os.path.exists(db_png_path):
                    os.link(db_png_path, out_dir + "/pathways/" + path + ".png")
                cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                      + db_png_path + ' ' + out_dir + "/pathways/" + path + ".pdf"
                cmd1_list.append(cmd)

        with open(self.work_dir + "/cmd0.list", "w") as fw:
            for i in range(len(cmd0_list)):
                fw.write(cmd0_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd0.list")
        cmd0_obj = self.add_command("cmd0", cmd, ignore_error=True).run()
        self.wait(cmd0_obj)

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
                    try:
                        color = link.split("%09")[1]  # tomato
                    except:
                        continue
                    dict[ko] = color
                self.map_dict[pathway] = dict  # self.map_dict["map05340"] = {"K10887": "yellow"}

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
                return "blue"
            else:
                return False
        elif len(self.category) == 2:
            # lst = list(self.category.keys())  # 基因集列表
            lst = self.proteinset_list
            # lst.sort()
            if ko in self.category[lst[0]] and ko in self.category[lst[1]]:
                return "pink"
            elif ko in self.category[lst[0]]:
                return "blue"
            elif ko in self.category[lst[1]]:
                return "red"
            else:
                return False
