# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from collections import defaultdict
import gridfs
import re
import tarfile
from bson import ObjectId
import time
from mbio.packages.ref_rna_v2.kegg_regulate import KeggRegulate
from biocluster.config import Config
import subprocess
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
import shutil
import unittest

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
            {"name": "geneset_kegg", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "task_id", "type": "string","default": ""},
            {"name": "kegg_table", "type": "infile", "format": "ref_rna_v2.kegg_table"},
            {"name": "kegg_table2", "type": "string", "default": ""},
            {"name": "geneset_id", "type": "string"},
            {"name": "annot_result", "type": "string"},
            {"name": "level", "type": "string"},
            {"name": "background_links", "type": "string", "default": ""},  # 底图的地址信息，add_info
            {"name": "type", "type": "string", "default": "origin"},  # 取最新的注释表还是原来的注释表
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
        if self.option('kegg_version') == "201909":
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Annotation/other2019/kegg201909/html')
        else:
            self.db_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.species_abr = 'map'
        self.gene2set = dict()
        self.geneset_ko = list()
        self.geneset_list = list()
        self.genome_path = self.config.SOFTWARE_DIR + "/database/KEGG/genome2.xls"
        self.map_dict = {}
        if self.option("task_id"):
            self.geneset_id_number = len(self.option('geneset_id').split(','))
        else:
            self.geneset_id_number = 1
            # self.colors = ['#DC143C', '#008000'] if self.option('source') == 'diff_exp' else ["#0000cd", "#ff8000"]
        # self.colors = ['#FF0000', '#00FF00']
        self.colors = ['#FF0000', '#0000FF']

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
        self.generate_ko_txt_dir()
        self.generate_new_pics()
        self.packing_press()
        self.end()

    def packing_press(self):
        if os.path.exists(os.path.join(self.output_dir,"pathways.tar.gz")):
            os.remove(os.path.join(self.output_dir,"pathways.tar.gz"))
        cmd = "tar -zcvf {} -C {} . ".format(os.path.join(self.output_dir,"pathways.tar.gz"), os.path.join(self.output_dir,"pathways"))
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('kegg图片打包压缩完成')
        except subprocess.CalledProcessError:
            self.logger.info('kegg图片打包压缩失败')
            self.set_error("kegg图片打包压缩失败失败")

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
        geneset_id = self.option("geneset_id")
        self.mongo_db = Config().get_mongo_client(mtype="whole_transcriptome")[Config().get_mongo_dbname("whole_transcriptome")]
        fs = gridfs.GridFS(self.mongo_db)
        annotation_collection = self.mongo_db["annotation_kegg"]
        geneset_collection = self.mongo_db["geneset"]
        if "," in geneset_id:
            geneset_id = geneset_id.split(",")[0]
        result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not result:
            self.set_error('geneset with main_id: %s was not found', variables = (geneset_id), code =  "33706401")
        # task_id = result["task_id"]
        anno_type = result["level"]
        main_id = annotation_collection.find_one({"task_id": self.option('task_id')})["main_id"]
        kegg_level_collection = self.mongo_db["annotation_kegg_level"]
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

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            os.system('ln -s {} {}'.format(path.rstrip("/"), to_path.rstrip("/")))
        else:
            try:
                transfer = MultiFileTransfer()
                transfer.add_download(path, to_path)
                transfer.perform()
            except:
                self.set_error('file can not find %s', variables=(path), code="33706404")
        if not os.path.exists(to_path):
            self.set_error("文件下载失败")
        elif not os.listdir(to_path):
            self.set_error("下载的文件夹为空")
        else:
            return to_path

    def get_kegg_png_workflow(self):

        # 工作流跑的时候还需要判定一下物种
        import pandas as pd
        kegg_df = pd.read_csv(self.option('kegg_table').prop['path'], sep='\t').fillna('')
        paths = sum([p.split(';') for p in kegg_df['Paths'] if p], [])
        specie = paths[0]
        self.species_abr = filter(str.isalpha, specie)
        if self.species_abr.lower() == 'ko':
            self.species_abr = 'map'

        geneset_type = self.option("level")
        kegg_annot_output = self.option("annot_result")
        class_dir = 'allannot_class'
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
            pass
            # self.set_error('file can not find %s', variables=kegg_pathways, code="33706404")
        if os.path.exists("png"):
            shutil.rmtree("png")
        if not os.path.exists("raw_png"):
            self.set_error('raw_png下载失败')
        shutil.copytree("raw_png", "png")

    def get_kegg_png(self):
        geneset_id = self.option("geneset_id")
        self.mongo_db = Config().get_mongo_client(mtype="whole_transcriptome")[Config().get_mongo_dbname("whole_transcriptome")]
        kegg_collection = self.mongo_db["task"]
        geneset_collection = self.mongo_db["geneset"]
        if "," in geneset_id:
            geneset_id = geneset_id.split(",")[0]
        result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not result:
            self.set_error('geneset with main_id: %s was not found', variables = (geneset_id), code =  "33706402")
        # task_id = result["task_id"]
        geneset_type = result["level"]  #全转录组字段名由type改成了level

        # kegg_annot_output = kegg_collection.find_one({"task_id": self.option('task_id'),'status': "end"})['sub_output']["long"]
        kegg_annot_output = kegg_collection.find_one({"task_id": self.option('task_id')})['sub_output']["long"]
        # class_dir = 'allannot_class'
        # if kegg_collection.find_one({"task_id": self.option('task_id'), 'status': "end"})['has_new'] == False:
        #     class_dir = 'refannot_class'

        # if geneset_type == 'G':
        #     kegg_pathways = os.path.join(kegg_annot_output,'annotation/allannot_class/kegg','kegg_pathway_gene_dir')
        # else:
        #     kegg_pathways = os.path.join(kegg_annot_output,'annotation/allannot_class/kegg', 'kegg_pathway_tran_dir')
        # anno_path = self.work_dir + "/png/"
        #
        # if not kegg_pathways.endswith("/"):
        #     kegg_pathways += "/"
        # anno_path = self.work_dir + "/raw_png/"
        # self.download_s3_file(kegg_pathways, anno_path)
        # shutil.copytree("raw_png", "png")

        if geneset_type == 'G':
            kegg_pathways = os.path.join(kegg_annot_output, 'annotation/allannot_class/kegg', 'kegg_pathway_gene_dir')
        else:
            kegg_pathways = os.path.join(kegg_annot_output, 'annotation/allannot_class/kegg', 'kegg_pathway_tran_dir')
        if exists(kegg_pathways):
            anno_path = self.work_dir + "/png/"

            if not kegg_pathways.endswith("/"):
                kegg_pathways += "/"
            anno_path = self.work_dir + "/raw_png/"
            self.download_s3_file(kegg_pathways, anno_path)
            shutil.copytree("raw_png", "png")

        else:
            kegg_pathways = os.path.join(kegg_annot_output, "annotation.tar.gz")
            rawfile = self.work_dir + "/annotation.tar.gz"
            download(kegg_pathways, rawfile)
            with tarfile.open(rawfile, 'r:gz') as tar:
                def is_within_directory(directory, target):
                    
                    abs_directory = os.path.abspath(directory)
                    abs_target = os.path.abspath(target)
                
                    prefix = os.path.commonprefix([abs_directory, abs_target])
                    
                    return prefix == abs_directory
                
                def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                
                    for member in tar.getmembers():
                        member_path = os.path.join(path, member.name)
                        if not is_within_directory(path, member_path):
                            raise Exception("Attempted Path Traversal in Tar File")
                
                    tar.extractall(path, members, numeric_owner=numeric_owner) 
                    
                
                safe_extract(tar, self.work_dir)
            if geneset_type == 'G':
                anno_path = os.path.join(self.work_dir, 'annotation/allannot_class/kegg', 'kegg_pathway_gene_dir')
            else:
                anno_path = os.path.join(self.work_dir, 'annotation/allannot_class/kegg', 'kegg_pathway_tran_dir')
            shutil.copytree(anno_path, "png")

    def get_dicts(self):
        gene_kegg_gene = dict()
        if self.species_abr == 'map':
            ko_genes, path_ko = self.option('kegg_table').get_pathway_koid(self.option('kegg_table2'), kegg_version=self.option('kegg_version'))
        else:
            #区分比对到整个ko与单个物种的ko时有区别
            ko_genes, path_ko, gene_kegg_gene = self.option('kegg_table').get_pathway_koid2()
        geneset_ko = defaultdict(set)
        regulate_gene = {}
        gene2set = dict()
        geneset2_ko = list()
        if self.option('source') != 'diff_exp':
            with open(self.option("geneset_kegg"), "r") as f:
                for line in f:
                    line = line.strip().split("\t")
                    if len(line) == 2:
                        regulate_gene[line[0]] = line[1].split(",")
                        geneset_ko[line[0]] = []
                        self.geneset_list.append(line[0])
                        for key in ko_genes.keys():
                            for gene in regulate_gene[line[0]]:
                                if gene in ko_genes[key]:
                                    geneset_ko[line[0]].append(key)

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
                        geneset2_ko.append(set(geneset_ko[line[0]]))
        # 20191202添加，考虑到虽然来源于差异但是只有up或者down的情况
        else:
            with open(self.option("geneset_kegg"), "r") as f:
                for index, line in enumerate(f):
                    line = line.strip().split("\t")
                    if len(line) == 2:
                        regulate_gene[line[0]] = line[1].split(",")
                        geneset_ko[line[0]] = []
                        self.geneset_list.append(line[0])
                        for key in ko_genes.keys():
                            for gene in regulate_gene[line[0]]:
                                if gene in ko_genes[key]:
                                    geneset_ko[line[0]].append(key)

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
                        geneset2_ko.append(set(geneset_ko[line[0]]))
                    else:
                        if index == 0:
                            self.colors = ['#00FF00','#FF0000']
                        else:
                            self.colors = ['#FF0000', '#00FF00']

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
        KeggRegulate().get_regulate_table(ko_gene=self.ko_genes, path_ko=self.path_ko,
                                          regulate_gene=self.geneset_gene,
                                          output= self.work_dir + '/kegg_stat_tmp.xls',
                                          gene_list = self.geneset_list)
        self.logger.info("生成kegg_stat_tmp.xls文件完毕")

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
                        color_dict[gene].append(self.colors[0])  # 蓝色(未必)
                elif len(self.geneset_gene) ==2:
                    if gene in gene1_list and gene in gene2_list:
                        color_dict[gene].append(','.join(self.colors))  # 蓝色,大红(未必)
                    elif gene in gene1_list:
                        color_dict[gene].append(self.colors[0])  # 蓝色(未必)
                    elif gene in gene2_list:
                        color_dict[gene].append(self.colors[1])  # 大红(未必)
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
                ko = path.replace("map", "ko")
                cmd = "{} {} {} {} {} {} {}".format(self.r_path, self.map_path, path,
                                                 ko_path, out_dir + "/pathways/" + path + ".png",
                                                 self.db_path + path + ".kgml",
                                                 self.work_dir + "/png/" + path + ".png")
                if self.option('kegg_version'):
                    ve = self.option('kegg_version')
                else:
                    ve = "202003"
                map_html = KeggHtml(ve)
                map_html.color_fg = self.colors
                try:
                    map_html.run_gene_set(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html", self.gene2set)
                    if os.path.exists(self.work_dir + "/png/" + path + ".html.mark"):
                        map_html.run_gene_set_mark(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/png/" + path + ".html.mark", self.geneset_ko)
                    elif os.path.exists(self.work_dir + "/png/" + path + ".html"):
                        map_html.run_gene_set_mark_from_html(path + '.png', self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/png/" + path + ".html", self.geneset_ko)
                    else:
                        self.logger.info("{}注释结果中没有html 或 mark文件".format(path))
                    '''
                    # html文件里缺少信息暂时无法兼容
                    elif os.path.exists(self.work_dir + "/png/" + path + ".html"):
                        map_html.run_gene_set_mark_from_html(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/png/" + path + ".html", self.geneset_ko)
                    '''
                except Exception as e:
                    self.logger.info("kegg_html wrong {}".format(e))
                    self.logger.info("注释结果中没有kegghtml {}".format(path))
                cmd1_list.append(cmd)
                ## 去掉pnd转成pdf的步骤，运行太耗时
                # pdf_path = out_dir + "/pathways/" + path + ".pdf"
                # cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                #       + out_dir + "/pathways/" + path + ".png" + ' ' + pdf_path
                # cmd2_list.append(cmd)
            else:
                db_png_path =  self.work_dir + "/png/" + path + ".png"
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
                ko_link_all = []
                ko_link_bg = []
                if pathway in self.map_dict:  # 含有背景色
                    for ko in self.map_dict[pathway].keys():  # 对背景色中的所有项进行循环
                        # self.logger.info(ko)
                        if ko == "":
                            continue
                        if ko in ko_tmp:  # 基因ko显著富集
                            if self.get_color(ko):
                                ko_link_all.append(ko + "+{},{}%0d%0a".format(self.map_dict[pathway][ko], self.get_color(ko)))
                                # lnk += ko + "+{},{}%0d%0a".format(self.map_dict[pathway][ko], self.get_color(ko))  # 有背景色,前景色
                            else:
                                ko_link_bg.append(ko + "+{}%0d".format(self.map_dict[pathway][ko]))
                                # lnk += ko + "+{}%0d".format(self.map_dict[pathway][ko])  # 只有背景色
                        else:
                            ko_link_bg.append(ko + "+{}%0d".format(self.map_dict[pathway][ko]))
                            # lnk += ko + "+{}%0d".format(self.map_dict[pathway][ko])  # 只有背景色
                else:  # 只标边框
                    for ko in ko_tmp:
                        if ko == "":
                            continue
                        if self.get_color(ko):
                            ko_link_bg.append(ko + "+{},{}%0d%0a".format("white", self.get_color(ko)))
                            # lnk += ko + "+{},{}%0d%0a".format("white", self.get_color(ko))  # 只有背景色
                        else:
                            pass
                tmp[-1] = lnk + "".join(ko_link_all) + "".join(ko_link_bg)
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
            "name": "whole_transcriptome.geneset.kegg_class",
            "instant": False,
            "options": dict(
                geneset_kegg="/mnt/ilustre/users/sanger-dev/workspace/20191129/GenesetEnrich_tsg_36238_8733_6071/multi_geneset_list",
                #group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
                task_id="tsg_36238",
                kegg_table="/mnt/ilustre/users/sanger-dev/workspace/20191129/GenesetEnrich_tsg_36238_8733_6071/gene_kegg_table.xls",
                kegg_table2="/mnt/ilustre/users/sanger-dev/workspace/20191129/GenesetEnrich_tsg_36238_8733_6071/gene_kegg_level_table.xls",
                geneset_id="5de0859a17b2bf6c0d644c2a",
                background_links="/mnt/ilustre/users/sanger-dev/workspace/20191129/GenesetEnrich_tsg_36238_8733_6071/add_info.txt",
                # type="origin",
                #corr_main_id="test_0711",
                #scm="complete",
                #scd="correlation",
                #corr_method="pearson",
                #Draw_in_groups="no",
                #log_base="10",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
