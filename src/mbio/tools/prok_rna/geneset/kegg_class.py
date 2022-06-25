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
import tarfile
from mbio.packages.prok_rna.kegg_regulate import KeggRegulate
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.config import Config
from mbio.packages.rna.annot_config import AnnotConfig
import subprocess
import shutil

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
            {"name": "task_id", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "geneset_kegg", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "annot_result", "type": "string"},
            {"name": "kegg_table", "type": "infile", "format": "prok_rna.kegg_table"},
            {'name': 'kegg_version', 'type': 'string', 'default': "2017"},
            {"name": "geneset_id", "type": "string"},
            {"name": "background_links", "type": "string", "default": ""},  # 底图的地址信息，add_info
            {"name": "type", "type": "string", "default": "origin"},  # 取最新的注释表还是原来的注释表
            {"name": "source", "type": "string", "default": ""},
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
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))

        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map5.r"
        self.db_path = self.kegg_files_dict["html"] + "/"
        # if self.option('kegg_version') in ["201909", "202003"]:
        #     self.db_path = os.path.join(self.config.SOFTWARE_DIR,
        #                                 'database/Annotation/other{}/kegg{}/html/'.format(self.option('kegg_version')[0:4],
        #                                 self.option('kegg_version')))
        # else:
        #     self.db_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
        # self.db_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
        self.geneset_ko = list()

        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.map_dict = {}
        self.species_abr = 'map'
        self.gene2set = dict()
        self.genome_path = self.config.SOFTWARE_DIR + "/database/KEGG/genome2.xls"
        self.geneset_list = list()
        self.colors = ['#0000FF', '#FF0000']

    def run(self):
        """
        运行
        :return:
        """
        super(KeggClassTool, self).run()
        # self.get_kegg_pics()
        if self.option("task_id"):
            self.get_kegg_png()
        else:
            self.get_kegg_png_workflow()
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
                os.system('ln -s {} {}'.format(path.rstrip("/"), to_path.rstrip("/")))
        else:
            # try:
            #     transfer = MultiFileTransfer()
            #     transfer.add_download(path, to_path)
            #     transfer.perform()
            # except:
            #     self.set_error('file can not find {}'.format(path))
            transfer = MultiFileTransfer()
            transfer.add_download(path, to_path)
            transfer.perform()
        return to_path

    def get_kegg_pics(self):
        """该函数被get_kegg_png替换"""
        geneset_id = self.option("geneset_id")
        self.mongo_db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
        fs = gridfs.GridFS(self.mongo_db)
        annotation_collection = self.mongo_db["sg_annotation_kegg"]
        geneset_collection = self.mongo_db["sg_geneset"]
        if "," in geneset_id:
            geneset_id = geneset_id.split(",")[0]
        result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not result:
            self.set_error("geneset with main_id: %s was not found", variables = (geneset_id), code = "35003501")
        task_id = result["task_id"]
        anno_type = result["type"]
        main_id = annotation_collection.find_one({"task_id":self.option('task_id'), "type": self.option('type')})["main_id"]

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

    def get_kegg_png_workflow(self):

        kegg_annot_output = self.option("annot_result")
        class_dir = 'allannot_class'
        kegg_pathways = os.path.join(kegg_annot_output,  'kegg', 'pathways')
        anno_path = self.work_dir + "/png/"
        if not kegg_pathways.endswith("/"):
            kegg_pathways += "/"
        anno_path = self.work_dir + "/raw_png/"
        self.logger.info("lalalalalalal{}".format(kegg_pathways))
        try:
            self.download_s3_file(kegg_pathways, anno_path)
        except:
            pass
            # self.set_error('file can not find %s', variables=kegg_pathways, code="33706404")
        if os.path.exists("png"):
            shutil.rmtree("png")
        if not os.path.exists("raw_png"):
            self.set_error('raw_png下载失败')
        shutil.copytree("raw_png", "png")
        self.anno_path = self.work_dir + "/png"



    def get_kegg_png(self):
        geneset_id = self.option("geneset_id")
        self.mongo_db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
        kegg_collection = self.mongo_db["sg_annotation_stat"]
        geneset_collection = self.mongo_db["sg_geneset"]
        if "," in geneset_id:
            geneset_id = geneset_id.split(",")[0]
        result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not result:
            self.set_error("geneset with main_id: %s was not found", variables = (geneset_id), code = "35003502")
        task_id = result["task_id"]
        kegg_annot_output = kegg_collection.find_one({"task_id": self.option('task_id'), "type": self.option('type'), 'status': "end"})['result_dir']
        annot_params = kegg_collection.find_one({"task_id": self.option('task_id'), "type": self.option('type'), 'status': "end"})['params']
        params_dict = json.loads(annot_params)
        if params_dict.has_key('kegg_species'):
            self.get_genome_abr(params_dict['kegg_species'])

        kegg_pathways = os.path.join(kegg_annot_output, 'kegg')
        if not kegg_pathways.endswith("/"):
            kegg_pathways += "/"
        anno_path = self.work_dir + "/kegg/"
        self.download_s3_file(kegg_pathways, anno_path)
        if os.path.exists(anno_path + "pathways.tar.gz"):
            tar = tarfile.open(anno_path + "pathways.tar.gz")
            tar.extractall(path=self.work_dir + "/png")
            self.anno_path = self.work_dir + "/png"
        elif os.path.exists(anno_path):
            self.anno_path = anno_path
        else:
            self.set_error("无法获得kegg注释结果图片信息", code = "35003503")

    def get_dicts(self):
        gene_kegg_gene = dict()
        if self.species_abr == 'map':
            ko_genes, path_ko = self.option('kegg_table').get_pathway_koid(kegg_version=self.option('kegg_version'))
        else:
            ko_genes, path_ko, gene_kegg_gene = self.option('kegg_table').get_pathway_koid2()
        geneset_ko = defaultdict(set)
        regulate_gene = {}
        gene2set = dict()
        geneset2_ko = list()
        with open(self.option("geneset_kegg"), "r") as f:
            for line in f:
                line = line.strip("\n").split("\t")
                if line[1] != "":
                    regulate_gene[line[0]] = line[1].split(",")
                else:
                    regulate_gene[line[0]] = []
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
                                for kegg_gene in kegg_genes.split(","):
                                    geneset_ko[line[0]].append(kegg_gene)
                        else:
                            pass
        for i in gene2set.keys():
            self.gene2set[i] = ",".join(gene2set[i])
        self.logger.info("path is {}".format(";".join(path_ko.keys())))

        self.category = geneset_ko  # 基因集与ko的对应关系
        self.geneset_ko = geneset2_ko
        self.ko_genes= ko_genes  # ko与基因的对应关系
        self.path_ko = path_ko  # path小ko与基因大Ko的对应关系
        self.geneset_gene = regulate_gene  # 基因集与基因的对应关系

        self.gene_kegg_gene = gene_kegg_gene #基因与kegg基因的对应关系

    def get_kegg_stat_tmp(self):
        self.logger.info("生成kegg_stat_tmp.xls文件")
        if self.species_abr == 'map':
            KeggRegulate().get_regulate_table(ko_gene=self.ko_genes, path_ko=self.path_ko,
                                              regulate_gene=self.geneset_gene,
                                              output= self.work_dir + '/kegg_stat_tmp.xls')
        else:
            KeggRegulate().get_regulate_table2(ko_gene=self.ko_genes, path_ko=self.path_ko,
                                              regulate_gene=self.geneset_gene,
                                              gene_kegg_gene=self.gene_kegg_gene,
                                              output= self.work_dir + '/kegg_stat_tmp.xls')

        self.logger.info("生成kegg_stat_tmp.xls文件完毕")
    def generate_ko_txt_dir2(self):
        kegg_regulate=self.output_dir + '/kegg_stat.xls'
        out_dir=self.output_dir
        if not os.path.exists(out_dir + "/ko"):
            os.mkdir(out_dir + "/ko")
        f = open(kegg_regulate, "r")
        headline = f.readline()
        headers = headline.split("\t")
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
                # 表格基因集顺序可能顺序改变
                if headers[3].startswith(self.geneset_list[0]):
                    pass
                else:
                    gene1_list, gene2_list = gene2_list,gene1_list
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
        headline = f.readline()
        headers = headline.split("\t")
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
                # 表格基因集顺序可能顺序改变
                if headers[3].startswith(self.category.keys()[0]):
                    pass
                else:
                    gene1_list, gene2_list = gene2_list, gene1_list
            gene_list = []
            gene_list.extend(gene1_list)
            gene_list.extend(gene2_list)
            gene_list_unrepeat = list(set(gene_list))
            color_dict = {}
            for gene in gene_list_unrepeat:
                color_dict[gene] = []
                if len(self.geneset_gene) == 1:
                    if gene in gene1_list:
                        color_dict[gene].append("#0000cd")  # 蓝色
                elif len(self.geneset_gene) ==2:
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
        cmd1_list = []
        cmd2_list = []
        for path in path_list:
            ko_path = kos_path + "/" + path
            if os.path.exists(ko_path):
                # ko = path.replace(self.species_abr, "ko")
                cmd = "{} {} {} {} {} {} {}".format(self.r_path, self.map_path, path,
                                                 ko_path, out_dir + "/pathways/" + path + ".png",
                                                 self.db_path + path + ".kgml",
                                                 self.anno_path + "/" + path + ".png")
                map_html = KeggHtml(self.option('kegg_version'))
                '''
                try:
                    map_html.run_gene_set(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html", self.gene2set)
                except:
                    self.logger.info("注释结果中没有kegghtml {}".format(path))
                '''

                if os.path.exists(self.anno_path + "/" + path + ".html.mark"):
                    map_html.run_gene_set_mark(self.anno_path + "/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.anno_path + "/" + path + ".html.mark", self.geneset_ko)
                else:
                    self.logger.info("{}注释结果中没有html 或 mark文件".format(path))


                cmd1_list.append(cmd)
                pdf_path = out_dir + "/pathways/" + path + ".pdf"
                cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                      + out_dir + "/pathways/" + path + ".png" + ' ' + pdf_path
                cmd2_list.append(cmd)
            else:
                db_png_path =  self.work_dir + "/png/" + path + ".png"
                if os.path.exists(out_dir + "/pathways/" + path + ".png"):
                    os.remove(out_dir + "/pathways/" + path + ".png")
                if os.path.exists(db_png_path):
                    os.link(db_png_path, out_dir + "/pathways/" + path + ".png")
                    cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                          + db_png_path + ' ' + out_dir + "/pathways/" + path + ".pdf"
                    cmd1_list.append(cmd)
        self.logger.info("开始生成新kegg图片")
        with open(self.work_dir + "/cmd1.list", "w") as fw:
            for i in range(len(cmd1_list)):
                fw.write(cmd1_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd1.list")
        cmd1_obj = self.add_command("cmd1", cmd, ignore_error=True).run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("cmd1 list执行成功")
        '''
        with open(self.work_dir + "/cmd2.list", "w") as fw:
            for i in range(len(cmd2_list)):
                fw.write(cmd2_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd2.list")
        cmd2_obj = self.add_command("cmd2", cmd, ignore_error=True).run()
        self.wait(cmd2_obj)
        if cmd2_obj.return_code == 0:
            self.logger.info("cmd2 list执行成功")
        '''
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
            lst = self.geneset_list
            # lst = list(self.category.keys())  # 基因集列表
            # lst.sort()
            if ko in self.category[lst[0]] and ko in self.category[lst[1]]:
                return "pink"
            elif ko in self.category[lst[0]]:
                return "blue"
            elif ko in self.category[lst[1]]:
                return "red"
            else:
                return False
